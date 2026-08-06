"""Targeted read remapping and conservative provisional-junction support reporting."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Callable, Iterable
from pathlib import Path
from tempfile import TemporaryDirectory

from contigger.evidence.bam import BamEvidenceProvider
from contigger.exceptions import InputValidationError
from contigger.fasta import write_fasta_records
from contigger.models import JunctionRemappingRequest, TargetedJunctionEvidence
from contigger.utilities.subprocesses import CommandResult, find_executable, run_command

CommandRunner = Callable[[tuple[str, ...]], CommandResult]
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


class TargetedJunctionRemapper:
    """Use samtools-selected sample reads to test one supplied provisional junction."""

    def __init__(
        self,
        source: BamEvidenceProvider,
        *,
        minimap2: str | Path = "minimap2",
        preset: str = "sr",
        threads: int = 1,
        runner: CommandRunner = run_command,
    ) -> None:
        if preset not in {"sr", "map-ont", "map-hifi"}:
            raise InputValidationError(f"unsupported targeted-remapping preset: {preset}")
        if threads < 1:
            raise InputValidationError("targeted-remapping thread count must be at least one")
        executable = Path(minimap2)
        if executable.is_file():
            executable = executable.resolve()
        elif len(executable.parts) == 1:
            found = find_executable(str(minimap2))
            if found is None:
                raise InputValidationError(f"required executable is unavailable: {minimap2}")
            executable = found
        else:
            raise InputValidationError(f"required executable does not exist: {minimap2}")
        self.source = source
        self.executable = str(executable)
        self.preset = preset
        self.threads = threads
        self.runner = runner
        self._commands: list[tuple[str, ...]] = []

    @property
    def commands(self) -> tuple[tuple[str, ...], ...]:
        """Return exact minimap2 commands executed by this remapper."""
        return tuple(self._commands)

    def evaluate(self, request: JunctionRemappingRequest) -> TargetedJunctionEvidence:
        """Extract relevant reads, remap them, and report distinct junction spanners."""
        if request.sample != self.source.sample.sample:
            raise InputValidationError(
                f"junction request sample {request.sample!r} does not match provider sample "
                f"{self.source.sample.sample!r}"
            )
        self._commands.clear()
        source_command_count = len(self.source.commands)
        left = tuple(
            self.source.reads_near_end(
                request.left_contig_id, request.left_end, request.extraction_distance
            )
        )
        right = tuple(
            self.source.reads_near_end(
                request.right_contig_id, request.right_end, request.extraction_distance
            )
        )
        selected = tuple(sorted(set(left) | set(right)))
        evaluation_commands = list(self.source.commands[source_command_count:])
        minimap_version = self._run((self.executable, "--version")).stdout.strip()
        evaluation_commands.append(self.commands[-1])
        if not minimap_version:
            raise InputValidationError("minimap2 returned an empty version response")
        assert self.source.version is not None
        if not selected:
            return self._result(
                request,
                selected,
                (),
                (),
                minimap_version,
                tuple(evaluation_commands),
                ("no reads selected",),
            )
        with TemporaryDirectory(prefix="contigger-junction-") as directory:
            work = Path(directory)
            reference = work / "junction.fasta"
            reads = work / "reads.fastq"
            write_fasta_records((request.provisional_reference,), reference)
            source_command_count = len(self.source.commands)
            self.source.extract_reads(selected, reads)
            evaluation_commands.extend(self.source.commands[source_command_count:])
            command = (
                self.executable,
                "-a",
                "-x",
                self.preset,
                "-t",
                str(self.threads),
                str(reference),
                str(reads),
            )
            sam = self._run(command).stdout
            evaluation_commands.append(self.commands[-1])
        remapped, spanning = _score_sam(
            sam,
            request.provisional_reference.identifier,
            request.provisional_reference.length,
            request.junction_position,
            request.minimum_spanning_flank,
            selected,
        )
        diagnostics = (
            "junction-spanning alignments are evidence only and do not authorize a merge",
        )
        return self._result(
            request,
            selected,
            remapped,
            spanning,
            minimap_version,
            tuple(evaluation_commands),
            diagnostics,
        )

    def _run(self, arguments: tuple[str, ...]) -> CommandResult:
        self._commands.append(arguments)
        return self.runner(arguments)

    def _result(
        self,
        request: JunctionRemappingRequest,
        selected: tuple[str, ...],
        remapped: tuple[str, ...],
        spanning: tuple[str, ...],
        minimap_version: str,
        commands: tuple[tuple[str, ...], ...],
        diagnostics: tuple[str, ...],
    ) -> TargetedJunctionEvidence:
        assert self.source.version is not None
        return TargetedJunctionEvidence(
            sample=request.sample,
            left_contig_id=request.left_contig_id,
            right_contig_id=request.right_contig_id,
            provisional_reference_id=request.provisional_reference.identifier,
            provisional_reference_length=request.provisional_reference.length,
            provisional_reference_sha256=hashlib.sha256(
                request.provisional_reference.sequence.encode("ascii")
            ).hexdigest(),
            junction_position=request.junction_position,
            selected_read_names=selected,
            remapped_read_names=remapped,
            spanning_read_names=spanning,
            minimum_spanning_flank=request.minimum_spanning_flank,
            samtools_version=self.source.version,
            minimap2_version=minimap_version,
            commands=commands,
            diagnostics=diagnostics,
        )


def _score_sam(
    text: str,
    reference_id: str,
    reference_length: int,
    junction: int,
    flank: int,
    selected_names: Iterable[str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    selected = set(selected_names)
    remapped: set[str] = set()
    spanning: set[str] = set()
    observed_reference = False
    for line_number, line in enumerate(text.splitlines(), start=1):
        if line.startswith("@SQ\t"):
            header_fields = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
            if header_fields.get("SN") == reference_id:
                observed_reference = header_fields.get("LN") == str(reference_length)
            continue
        if not line or line.startswith("@"):
            continue
        record_fields = line.split("\t")
        if len(record_fields) < 11:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: expected at least 11 fields"
            )
        name, flag_text, target, position_text, _mapq, cigar = record_fields[:6]
        if name not in selected:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: unexpected read name {name!r}"
            )
        try:
            flag = int(flag_text)
            start = int(position_text) - 1
        except ValueError as error:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: invalid numeric field"
            ) from error
        if flag & 0x4 or flag & 0x900 or target == "*":
            continue
        if target != reference_id:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: unexpected reference {target!r}"
            )
        reference_span, aligned_blocks = _reference_alignment(cigar, line_number)
        end = start + reference_span
        if not 0 <= start < end <= reference_length:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: alignment is outside the provisional reference"
            )
        remapped.add(name)
        if any(
            start + block_start <= junction - flank and start + block_end >= junction + flank
            for block_start, block_end in aligned_blocks
        ):
            spanning.add(name)
    if not observed_reference:
        raise InputValidationError(
            "minimap2 SAM header does not contain the exact provisional reference name and length"
        )
    return tuple(sorted(remapped)), tuple(sorted(spanning))


def _reference_alignment(cigar: str, line_number: int) -> tuple[int, tuple[tuple[int, int], ...]]:
    if cigar == "*":
        return 0, ()
    parts = _CIGAR.findall(cigar)
    if not parts or "".join(length + operation for length, operation in parts) != cigar:
        raise InputValidationError(f"minimap2 SAM line {line_number}: malformed CIGAR {cigar!r}")
    position = 0
    block_start: int | None = None
    blocks: list[tuple[int, int]] = []
    for length_text, operation in parts:
        length = int(length_text)
        if length < 1:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: CIGAR operations must be positive"
            )
        if operation in "M=X":
            if block_start is None:
                block_start = position
            position += length
        elif operation in "DN":
            if block_start is not None:
                blocks.append((block_start, position))
                block_start = None
            position += length
        # Insertions and clipping consume no reference and do not bridge a reference gap.
    if block_start is not None:
        blocks.append((block_start, position))
    return position, tuple(blocks)
