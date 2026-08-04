"""Minimal minimap2 discovery, command construction, and PAF parsing scaffold."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

from contigger.exceptions import ExternalToolError, FeatureNotImplementedError, InputValidationError
from contigger.models import AlignmentHit, Orientation, SequenceRecord
from contigger.utilities.subprocesses import find_executable, run_command


class Minimap2Aligner:
    """A replaceable minimap2 adapter; complete workflow is intentionally deferred."""

    def __init__(self, executable: Path | None = None, threads: int = 1) -> None:
        self.executable = executable or find_executable("minimap2")
        self.threads = threads
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
            "asm20",
            "-t",
            str(self.threads),
            str(target),
            str(query),
        )


def parse_paf_line(line: str) -> AlignmentHit:
    """Parse the required twelve PAF fields and selected optional integer tags."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 12:
        raise InputValidationError("PAF record requires at least 12 tab-separated fields")
    try:
        orientation = Orientation(fields[4])
        tags = {
            field[:2]: field[5:] for field in fields[12:] if len(field) >= 6 and field[2:5] == ":i:"
        }
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
            alignment_score=int(tags["AS"]) if "AS" in tags else None,
            supporting_seeds=int(tags["cm"]) if "cm" in tags else None,
        )
    except (ValueError, KeyError) as error:
        raise InputValidationError(f"invalid PAF record: {error}") from error
