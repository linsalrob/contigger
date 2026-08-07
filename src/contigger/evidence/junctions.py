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
from contigger.textio import open_text
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
            request.minimum_mapping_quality,
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
            technology=self.source.sample.technology or "unknown",
            remapping_preset=self.preset,
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
            minimum_mapping_quality=request.minimum_mapping_quality,
            samtools_version=self.source.version,
            minimap2_version=minimap_version,
            commands=commands,
            diagnostics=diagnostics,
        )


class FastqJunctionRemapper:
    """Remap one sample's explicit FASTQ to a supplied provisional junction."""

    def __init__(
        self,
        *,
        minimap2: str | Path = "minimap2",
        preset: str = "map-ont",
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
        self.executable = str(executable)
        self.preset = preset
        self.threads = threads
        self.runner = runner
        self._version: str | None = None

    def evaluate(
        self,
        request: JunctionRemappingRequest,
        *,
        sample: str,
        technology: str,
        fastq: Path,
    ) -> TargetedJunctionEvidence:
        """Report distinct remapped and continuously junction-spanning FASTQ reads."""
        if request.sample != sample:
            raise InputValidationError(
                f"junction request sample {request.sample!r} does not match FASTQ sample {sample!r}"
            )
        selected = _read_fastq_names(fastq)
        version_command = (self.executable, "--version")
        commands: list[tuple[str, ...]] = []
        if self._version is None:
            self._version = self.runner(version_command).stdout.strip()
            commands.append(version_command)
            if not self._version:
                raise InputValidationError("minimap2 returned an empty version response")
        with TemporaryDirectory(prefix="contigger-fastq-junction-") as directory:
            reference = Path(directory) / "junction.fasta"
            write_fasta_records((request.provisional_reference,), reference)
            command = (
                self.executable,
                "-c",
                "--secondary=no",
                "-x",
                self.preset,
                "-t",
                str(self.threads),
                str(reference),
                str(fastq),
            )
            paf = self.runner(command).stdout
            commands.append(command)
        remapped, spanning = _score_paf_junction(
            paf,
            request.provisional_reference.identifier,
            request.provisional_reference.length,
            request.junction_position,
            request.minimum_spanning_flank,
            request.minimum_mapping_quality,
            selected,
        )
        return TargetedJunctionEvidence(
            sample=sample,
            technology=technology,
            remapping_preset=self.preset,
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
            minimum_mapping_quality=request.minimum_mapping_quality,
            samtools_version="not used",
            minimap2_version=self._version,
            commands=tuple(commands),
            diagnostics=(
                "checked-in sample FASTQ remapping is benchmark evidence only and does not "
                "authorize a merge",
            ),
        )


def _read_fastq_names(path: Path) -> tuple[str, ...]:
    names: list[str] = []
    seen: set[str] = set()
    with open_text(path, encoding="utf-8", newline="") as handle:
        line_number = 0
        while header := handle.readline():
            line_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise InputValidationError(f"{path}:{line_number}: incomplete FASTQ record")
            if not header.startswith("@"):
                raise InputValidationError(
                    f"{path}:{line_number}: FASTQ header must start with '@'"
                )
            if not plus.startswith("+"):
                raise InputValidationError(
                    f"{path}:{line_number + 2}: FASTQ separator must start with '+'"
                )
            name = header[1:].strip().split(maxsplit=1)[0]
            if not name:
                raise InputValidationError(f"{path}:{line_number}: empty FASTQ read name")
            if name in seen:
                raise InputValidationError(
                    f"{path}:{line_number}: duplicate FASTQ read name {name!r}"
                )
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise InputValidationError(
                    f"{path}:{line_number}: FASTQ sequence and quality lengths differ"
                )
            seen.add(name)
            names.append(name)
            line_number += 3
    if not names:
        raise InputValidationError(f"{path}: FASTQ file is empty")
    return tuple(sorted(names))


def _score_paf_junction(
    text: str,
    reference_id: str,
    reference_length: int,
    junction: int,
    flank: int,
    minimum_mapping_quality: int,
    selected_names: Iterable[str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    selected = set(selected_names)
    remapped: set[str] = set()
    spanning: set[str] = set()
    for line_number, line in enumerate(text.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) < 12:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: expected at least 12 fields"
            )
        name, target = fields[0], fields[5]
        if name not in selected:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: unexpected read name {name!r}"
            )
        if target != reference_id:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: unexpected reference {target!r}"
            )
        try:
            target_length = int(fields[6])
            start = int(fields[7])
            end = int(fields[8])
            mapping_quality = int(fields[11])
        except ValueError as error:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: invalid numeric field"
            ) from error
        if not 0 <= mapping_quality <= 255:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: mapping quality must be between 0 and 255"
            )
        if target_length != reference_length or not 0 <= start < end <= reference_length:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: alignment is outside the provisional reference"
            )
        alignment_types = [field[5:] for field in fields[12:] if field.startswith("tp:A:")]
        if len(alignment_types) != 1:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: expected exactly one tp:A alignment-type tag"
            )
        if alignment_types[0] != "P":
            continue
        if minimum_mapping_quality and (
            mapping_quality == 255 or mapping_quality < minimum_mapping_quality
        ):
            continue
        cigars = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
        if len(cigars) != 1:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: expected exactly one cg:Z CIGAR tag"
            )
        reference_span, blocks = _reference_alignment(cigars[0], line_number)
        if reference_span != end - start:
            raise InputValidationError(
                f"minimap2 PAF line {line_number}: CIGAR reference span does not match coordinates"
            )
        remapped.add(name)
        if any(
            start + block_start <= junction - flank and start + block_end >= junction + flank
            for block_start, block_end in blocks
        ):
            spanning.add(name)
    return tuple(sorted(remapped)), tuple(sorted(spanning))


def _score_sam(
    text: str,
    reference_id: str,
    reference_length: int,
    junction: int,
    flank: int,
    minimum_mapping_quality: int,
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
        name, flag_text, target, position_text, mapq_text, cigar = record_fields[:6]
        if name not in selected:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: unexpected read name {name!r}"
            )
        try:
            flag = int(flag_text)
            start = int(position_text) - 1
            mapping_quality = int(mapq_text)
        except ValueError as error:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: invalid numeric field"
            ) from error
        if flag & 0x4 or flag & 0x900 or target == "*":
            continue
        if not 0 <= mapping_quality <= 255:
            raise InputValidationError(
                f"minimap2 SAM line {line_number}: mapping quality must be between 0 and 255"
            )
        if minimum_mapping_quality and (
            mapping_quality == 255 or mapping_quality < minimum_mapping_quality
        ):
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
