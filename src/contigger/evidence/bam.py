"""Sample-scoped samtools BAM/CRAM validation and source pileups."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path

from contigger.exceptions import FeatureNotImplementedError, InputValidationError
from contigger.fasta import read_fasta
from contigger.models import BaseEvidence, JoinEvidence, SampleInput
from contigger.utilities.subprocesses import CommandResult, find_executable, run_command

CommandRunner = Callable[[tuple[str, ...]], CommandResult]


class BamEvidenceProvider:
    """Samtools-backed evidence provider scoped to exactly one source sample."""

    def __init__(
        self,
        sample: SampleInput,
        *,
        executable: str | Path = "samtools",
        runner: CommandRunner = run_command,
    ) -> None:
        if sample.bam is None:
            raise InputValidationError(f"sample {sample.sample!r} has no BAM/CRAM input")
        if not sample.bam.is_file():
            raise InputValidationError(f"BAM/CRAM does not exist: {sample.bam}")
        if sample.bam.suffix.lower() not in {".bam", ".cram"}:
            raise InputValidationError(f"alignment input must end in .bam or .cram: {sample.bam}")
        resolved = Path(executable)
        if resolved.is_file():
            resolved = resolved.resolve()
        elif len(resolved.parts) == 1:
            found = find_executable(str(executable))
            if found is None:
                raise InputValidationError(f"required executable is unavailable: {executable}")
            resolved = found
        else:
            raise InputValidationError(f"required executable does not exist: {executable}")
        self._sample = sample
        self._executable = str(resolved)
        self._runner = runner
        self._references: dict[str, int] | None = None
        self._version: str | None = None
        self._commands: list[tuple[str, ...]] = []

    @property
    def sample(self) -> SampleInput:
        """Return the provider's sample input."""
        return self._sample

    @property
    def version(self) -> str | None:
        """Return the captured samtools version after validation."""
        return self._version

    @property
    def commands(self) -> tuple[tuple[str, ...], ...]:
        """Return every exact samtools argument array executed by this provider."""
        return tuple(self._commands)

    def validate_source(self) -> tuple[str, ...]:
        """Validate index, integrity, and exact reference names/lengths against FASTA."""
        bam = self._require_bam()
        index = _alignment_index(bam)
        if index is None:
            raise InputValidationError(f"BAM/CRAM lacks an adjacent index: {bam}")
        version = self._run((self._executable, "--version")).stdout.splitlines()
        if not version:
            raise InputValidationError("samtools returned an empty version response")
        self._version = version[0].strip()
        self._run((self._executable, "quickcheck", "-v", str(bam)))
        self._run((self._executable, "idxstats", str(bam)))
        view_arguments = [self._executable, "view", "-H"]
        if bam.suffix.lower() == ".cram":
            view_arguments.extend(("-T", str(self.sample.contigs)))
        view_arguments.append(str(bam))
        header = self._run(tuple(view_arguments)).stdout
        observed = _parse_sq_header(header, bam)
        expected = {
            record.original_identifier: record.length
            for record in read_fasta(self.sample.contigs, self.sample.sample)
        }
        if observed != expected:
            missing = sorted(set(expected) - set(observed))
            extra = sorted(set(observed) - set(expected))
            if missing:
                detail = f"missing reference {missing[0]!r}"
            elif extra:
                detail = f"unexpected reference {extra[0]!r}"
            else:
                mismatch = next(
                    name for name in sorted(expected) if expected[name] != observed[name]
                )
                detail = (
                    f"reference {mismatch!r} length {observed[mismatch]} does not match "
                    f"FASTA length {expected[mismatch]}"
                )
            raise InputValidationError(
                f"sample {self.sample.sample!r} BAM/CRAM references do not match FASTA: {detail}"
            )
        self._references = observed
        return tuple(sorted(observed))

    def pileup(self, contig_id: str, start: int, end: int) -> Iterable[BaseEvidence]:
        """Return deterministic evidence for a zero-based, half-open source interval."""
        if self._references is None:
            self.validate_source()
        assert self._references is not None
        if contig_id not in self._references:
            raise InputValidationError(f"unknown BAM/CRAM source reference: {contig_id}")
        if not 0 <= start <= end <= self._references[contig_id]:
            raise InputValidationError(
                f"pileup interval is outside reference {contig_id!r}: [{start}, {end})"
            )
        if start == end:
            return ()
        region = f"{contig_id}:{start + 1}-{end}"
        arguments = [self._executable, "mpileup", "-aa", "-s"]
        if self._require_bam().suffix.lower() == ".cram":
            arguments.extend(("-f", str(self.sample.contigs)))
        arguments.extend(("-r", region, str(self._require_bam())))
        result = self._run(tuple(arguments))
        return tuple(_parse_pileup(result.stdout, self.sample.sample, contig_id, start, end))

    def reads_near_end(self, contig_id: str, end: str, distance: int) -> Iterable[str]:
        """Reserve end-read extraction for the targeted-remapping milestone."""
        raise FeatureNotImplementedError("read extraction near contig ends is not implemented")

    def junction_evidence(self, left_contig_id: str, right_contig_id: str) -> JoinEvidence:
        """State explicitly that a source BAM cannot validate a new junction."""
        return JoinEvidence(
            sample=self.sample.sample,
            left_contig_id=left_contig_id,
            right_contig_id=right_contig_id,
            testable=False,
            diagnostics=(
                "existing source alignment cannot validate a newly constructed junction; "
                "targeted remapping is not implemented",
            ),
        )

    def _require_bam(self) -> Path:
        assert self.sample.bam is not None
        return self.sample.bam

    def _run(self, arguments: tuple[str, ...]) -> CommandResult:
        self._commands.append(arguments)
        return self._runner(arguments)


def _alignment_index(path: Path) -> Path | None:
    if path.suffix.lower() == ".bam":
        candidates = [Path(f"{path}.bai"), path.with_suffix(".bai")]
    else:
        candidates = [Path(f"{path}.crai"), path.with_suffix(".crai")]
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def _parse_sq_header(header: str, path: Path) -> dict[str, int]:
    references: dict[str, int] = {}
    for line_number, line in enumerate(header.splitlines(), start=1):
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        name = fields.get("SN")
        length_text = fields.get("LN")
        if not name or not length_text:
            raise InputValidationError(f"{path}: header line {line_number}: @SQ requires SN and LN")
        if name in references:
            raise InputValidationError(f"{path}: header line {line_number}: duplicate @SQ {name!r}")
        try:
            length = int(length_text)
        except ValueError as error:
            raise InputValidationError(
                f"{path}: header line {line_number}: invalid @SQ length {length_text!r}"
            ) from error
        if length < 1:
            raise InputValidationError(
                f"{path}: header line {line_number}: @SQ length must be positive"
            )
        references[name] = length
    if not references:
        raise InputValidationError(f"{path}: BAM/CRAM header contains no @SQ references")
    return references


def _parse_pileup(
    text: str, sample: str, contig_id: str, start: int, end: int
) -> Iterable[BaseEvidence]:
    seen: set[int] = set()
    for line_number, line in enumerate(text.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != 7:
            raise InputValidationError(f"samtools mpileup line {line_number}: expected 7 fields")
        reference, position_text, reference_base, depth_text, bases, qualities, mappings = fields
        try:
            position = int(position_text) - 1
            depth = int(depth_text)
        except ValueError as error:
            raise InputValidationError(
                f"samtools mpileup line {line_number}: invalid numeric field"
            ) from error
        if depth < 0:
            raise InputValidationError(
                f"samtools mpileup line {line_number}: depth cannot be negative"
            )
        if reference != contig_id or not start <= position < end or position in seen:
            raise InputValidationError(
                f"samtools mpileup line {line_number}: unexpected reference or position"
            )
        seen.add(position)
        if depth == 0:
            alleles: dict[str, int] = {}
            base_scores: list[int] = []
            mapping_scores: list[int] = []
        else:
            alleles = _pileup_alleles(bases, reference_base, line_number)
            base_scores = [ord(value) - 33 for value in qualities]
            mapping_scores = [ord(value) - 33 for value in mappings]
        if len(base_scores) != len(mapping_scores):
            raise InputValidationError(
                f"samtools mpileup line {line_number}: base and mapping quality counts differ"
            )
        if any(score < 0 for score in base_scores + mapping_scores):
            raise InputValidationError(
                f"samtools mpileup line {line_number}: invalid quality encoding"
            )
        yield BaseEvidence(
            sample,
            contig_id,
            position,
            dict(sorted(alleles.items())),
            depth,
            sum(base_scores) / len(base_scores) if base_scores else None,
            sum(mapping_scores) / len(mapping_scores) if mapping_scores else None,
        )


def _pileup_alleles(bases: str, reference: str, line_number: int) -> dict[str, int]:
    counts: dict[str, int] = {}
    index = 0
    while index < len(bases):
        symbol = bases[index]
        if symbol == "^":
            if index + 1 >= len(bases):
                raise InputValidationError(
                    f"samtools mpileup line {line_number}: truncated read start"
                )
            index += 2
            continue
        if symbol == "$":
            index += 1
            continue
        if symbol in "+-":
            length_start = index + 1
            length_end = length_start
            while length_end < len(bases) and bases[length_end].isdigit():
                length_end += 1
            if length_end == length_start:
                raise InputValidationError(f"samtools mpileup line {line_number}: invalid indel")
            indel_length = int(bases[length_start:length_end])
            index = length_end + indel_length
            if index > len(bases):
                raise InputValidationError(f"samtools mpileup line {line_number}: truncated indel")
            continue
        if symbol in ".,":
            allele = reference.upper()
        elif symbol.upper() in {"A", "C", "G", "T", "N", "*"}:
            allele = symbol.upper()
        elif symbol in "<>":
            index += 1
            continue
        else:
            raise InputValidationError(
                f"samtools mpileup line {line_number}: unsupported base symbol {symbol!r}"
            )
        counts[allele] = counts.get(allele, 0) + 1
        index += 1
    return counts
