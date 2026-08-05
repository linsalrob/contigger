"""Safe minimap2 execution and strict streaming PAF parsing."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path

from contigger.exceptions import ExternalToolError, FeatureNotImplementedError, InputValidationError
from contigger.models import AlignmentHit, AlignmentType, Orientation, SequenceRecord
from contigger.utilities.subprocesses import find_executable, run_command


class Minimap2Aligner:
    """A replaceable path-based minimap2 assembly alignment adapter."""

    def __init__(
        self, executable: Path | None = None, threads: int = 1, preset: str = "asm20"
    ) -> None:
        if threads < 1:
            raise InputValidationError("minimap2 thread count must be at least one")
        if preset not in {"asm5", "asm10", "asm20"}:
            raise InputValidationError(f"unsupported minimap2 assembly preset: {preset}")
        self.executable = executable or find_executable("minimap2")
        self.threads = threads
        self.preset = preset
        self._last_command: tuple[str, ...] | None = None
        self._version: str | None = None

    @property
    def tool_name(self) -> str:
        """Return the adapter tool name."""
        return "minimap2"

    @property
    def tool_version(self) -> str:
        """Detect and cache the minimap2 version."""
        if self.executable is None:
            raise ExternalToolError("minimap2 executable not found")
        if self._version is None:
            command = (str(self.executable), "--version")
            self._last_command = command
            self._version = run_command(command).stdout.strip()
        return self._version

    @property
    def last_command(self) -> tuple[str, ...] | None:
        """Return the exact most recently executed command."""
        return self._last_command

    def build_index(self, targets: Sequence[SequenceRecord], index_path: Path) -> Path:
        """Reserve the indexing interface without reporting false completion."""
        raise FeatureNotImplementedError("minimap2 index construction is not implemented")

    def align(
        self, queries: Iterable[SequenceRecord], targets: Iterable[SequenceRecord]
    ) -> Iterable[AlignmentHit]:
        """Reserve the alignment interface without reporting biological results."""
        raise FeatureNotImplementedError("minimap2 alignment workflow is not implemented")

    def command_for_paths(self, target: Path, query: Path) -> tuple[str, ...]:
        """Construct a safe, provenance-ready assembly alignment command."""
        if self.executable is None:
            raise ExternalToolError("minimap2 executable not found")
        return (
            str(self.executable),
            "-x",
            self.preset,
            "-t",
            str(self.threads),
            str(target),
            str(query),
        )

    def align_paths(self, target: Path, query: Path) -> Iterator[AlignmentHit]:
        """Run minimap2 for two FASTA paths and parse its captured PAF output."""
        if not target.is_file():
            raise InputValidationError(f"target FASTA does not exist: {target}")
        if not query.is_file():
            raise InputValidationError(f"query FASTA does not exist: {query}")
        # Detect the version before alignment so both facts are available for provenance.
        _ = self.tool_version
        command = self.command_for_paths(target, query)
        self._last_command = command
        result = run_command(command)
        yield from parse_paf(result.stdout.splitlines())


def parse_paf_line(line: str) -> AlignmentHit:
    """Parse one PAF record, validating required fields and supported tags."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 12:
        raise InputValidationError("PAF record requires at least 12 tab-separated fields")
    if not fields[0] or not fields[5]:
        raise InputValidationError("PAF query and target names cannot be empty")
    try:
        orientation = Orientation(fields[4])
        tags = _parse_optional_tags(fields[12:])
        return AlignmentHit(
            query_id=fields[0],
            query_length=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            orientation=orientation,
            target_id=fields[5],
            target_length=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            matching_bases=int(fields[9]),
            alignment_block_length=int(fields[10]),
            mapping_quality=int(fields[11]),
            alignment_score=_integer_tag(tags, "AS"),
            supporting_seeds=_integer_tag(tags, "cm"),
            chaining_score=_integer_tag(tags, "s1"),
            secondary_chaining_score=_integer_tag(tags, "s2"),
            alignment_type=_alignment_type_tag(tags),
        )
    except (ValueError, KeyError, InputValidationError) as error:
        raise InputValidationError(f"invalid PAF record: {error}") from error


def parse_paf(lines: Iterable[str]) -> Iterator[AlignmentHit]:
    """Yield PAF hits from text lines, ignoring only blank lines.

    Errors include the one-based physical source line, including blank lines.
    """
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        try:
            yield parse_paf_line(line)
        except InputValidationError as error:
            raise InputValidationError(f"PAF line {line_number}: {error}") from error


def _parse_optional_tags(fields: Sequence[str]) -> dict[str, tuple[str, str]]:
    tags: dict[str, tuple[str, str]] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) != 3 or len(parts[0]) != 2 or len(parts[1]) != 1:
            raise InputValidationError(f"malformed optional tag {field!r}")
        name, tag_type, value = parts
        if name in tags:
            raise InputValidationError(f"duplicate optional tag {name!r}")
        tags[name] = (tag_type, value)
    return tags


def _integer_tag(tags: dict[str, tuple[str, str]], name: str) -> int | None:
    if name not in tags:
        return None
    tag_type, value = tags[name]
    if tag_type != "i":
        raise InputValidationError(f"tag {name!r} must have type 'i'")
    return int(value)


def _alignment_type_tag(tags: dict[str, tuple[str, str]]) -> AlignmentType | None:
    if "tp" not in tags:
        return None
    tag_type, value = tags["tp"]
    if tag_type != "A":
        raise InputValidationError("tag 'tp' must have type 'A'")
    return AlignmentType(value)
