"""Safe minimap2 execution and strict streaming PAF parsing."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

from contigger.exceptions import ExternalToolError, InputValidationError
from contigger.fasta import write_fasta_records
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
        """Build or safely reuse a content-validated minimap2 target index."""
        target_records = tuple(sorted(targets, key=lambda record: record.identifier))
        if not target_records:
            raise InputValidationError("minimap2 index requires at least one target")
        _validate_unique_identifiers(target_records, "index targets")
        metadata_path = _index_metadata_path(index_path)
        if index_path.exists() or metadata_path.exists():
            if not index_path.is_file() or not metadata_path.is_file():
                raise InputValidationError(f"incomplete minimap2 index or metadata at {index_path}")
            expected = self._index_metadata(target_records)
            try:
                observed = json.loads(metadata_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                raise InputValidationError(
                    f"cannot read minimap2 index metadata {metadata_path}: {error}"
                ) from error
            if observed != expected:
                raise InputValidationError(
                    f"minimap2 index metadata does not match requested targets: {index_path}"
                )
            return index_path
        expected = self._index_metadata(target_records)
        if self.executable is None:
            raise ExternalToolError("minimap2 executable not found")
        try:
            index_path.parent.mkdir(parents=True, exist_ok=True)
            with TemporaryDirectory(prefix=".contigger-index-", dir=index_path.parent) as directory:
                target_path = Path(directory) / "targets.fasta"
                temporary_index = Path(directory) / "target.mmi"
                write_fasta_records(target_records, target_path)
                command = (
                    str(self.executable),
                    "-x",
                    self.preset,
                    "-d",
                    str(temporary_index),
                    str(target_path),
                )
                self._last_command = command
                run_command(command)
                if not temporary_index.is_file():
                    raise ExternalToolError(
                        f"minimap2 did not create requested index: {index_path}"
                    )
                temporary_index.replace(index_path)
            metadata_path.write_text(
                json.dumps(expected, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
        except OSError as error:
            raise InputValidationError(
                f"cannot create minimap2 index {index_path}: {error}"
            ) from error
        return index_path

    def align_indexed(
        self,
        queries: Iterable[SequenceRecord],
        targets: Sequence[SequenceRecord],
        index_path: Path,
    ) -> Iterable[AlignmentHit]:
        """Align a query batch to an index proven to contain the supplied targets."""
        query_records = tuple(queries)
        target_records = tuple(targets)
        if not query_records:
            raise InputValidationError("indexed minimap2 alignment requires queries")
        _validate_unique_identifiers(query_records, "indexed queries")
        self.build_index(target_records, index_path)
        allowed_targets = {record.identifier for record in target_records}
        with TemporaryDirectory(prefix="contigger-align-") as directory:
            query_path = Path(directory) / "queries.fasta"
            write_fasta_records(query_records, query_path)
            for hit in self.align_paths(index_path, query_path):
                if hit.target_id not in allowed_targets:
                    raise InputValidationError(
                        f"minimap2 index returned an unexpected target identifier: {hit.target_id}"
                    )
                yield hit

    def _index_metadata(self, targets: Sequence[SequenceRecord]) -> dict[str, object]:
        """Return deterministic index identity without embedding sequence content."""
        digest = hashlib.sha256()
        identities: list[dict[str, object]] = []
        for record in targets:
            digest.update(record.identifier.encode("utf-8"))
            digest.update(b"\0")
            digest.update(record.sequence.encode("ascii"))
            digest.update(b"\0")
            identities.append({"identifier": record.identifier, "length": record.length})
        return {
            "format": 1,
            "minimap2_version": self.tool_version,
            "preset": self.preset,
            "sequence_sha256": digest.hexdigest(),
            "targets": identities,
        }

    def align(
        self, queries: Iterable[SequenceRecord], targets: Iterable[SequenceRecord]
    ) -> Iterable[AlignmentHit]:
        """Align supplied typed records through safe temporary FASTA paths.

        Selective callers pass one query and one target. Multi-record inputs remain
        valid minimap2 batches but are never constructed by the selective executor.
        """
        query_records = tuple(queries)
        target_records = tuple(targets)
        if not query_records or not target_records:
            raise InputValidationError("minimap2 alignment requires queries and targets")
        with TemporaryDirectory(prefix="contigger-align-") as directory:
            temporary = Path(directory)
            query_path = temporary / "queries.fasta"
            target_path = temporary / "targets.fasta"
            write_fasta_records(query_records, query_path)
            write_fasta_records(target_records, target_path)
            yield from self.align_paths(target_path, query_path)

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


def _index_metadata_path(index_path: Path) -> Path:
    return index_path.with_name(index_path.name + ".json")


def _validate_unique_identifiers(records: Sequence[SequenceRecord], label: str) -> None:
    identifiers = [record.identifier for record in records]
    if len(set(identifiers)) != len(identifiers):
        raise InputValidationError(f"{label} require unique identifiers")


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
