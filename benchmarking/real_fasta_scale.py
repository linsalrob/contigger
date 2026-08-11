"""Build deterministic real-data FASTA scale fixtures without copying source data into Git.

The selected record sets are nested: a 10,000-record set is a subset of a
100,000-record set made with the same seed.  Selection uses a stable hash of the
FASTA identifier, so it is independent of source-record order and approximately
preserves the source contig-length distribution.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import json
import os
import platform
import sys
import time
from collections.abc import Iterator
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True, slots=True)
class FastaRecordLocation:
    """One FASTA record identity and its deterministic selection priority."""

    identifier: str
    priority: int


@dataclass(frozen=True, slots=True)
class SourceSummary:
    """Counts measured while scanning a source FASTA."""

    records: int
    bases: int


def _open_fasta(path: Path) -> TextIO:
    """Open plain or gzip-compressed FASTA text for streaming."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def _records(path: Path) -> Iterator[tuple[str, list[str]]]:
    """Yield raw FASTA header identifiers and sequence lines without materialising all data."""
    identifier: str | None = None
    lines: list[str] = []
    with _open_fasta(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, lines
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                identifier = header.split(maxsplit=1)[0]
                lines = [line]
            elif identifier is None:
                if line.strip():
                    raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
            else:
                lines.append(line)
    if identifier is None:
        raise ValueError(f"{path}: FASTA file is empty")
    yield identifier, lines


def _priority(identifier: str, seed: str) -> int:
    """Return a stable, seed-scoped 64-bit rank for a FASTA identifier."""
    payload = f"{seed}\0{identifier}".encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")


def select_identifiers(path: Path, count: int, seed: str) -> tuple[set[str], SourceSummary]:
    """Select exactly ``count`` lowest-rank identifiers and measure source size."""
    if count < 1:
        raise ValueError("record count must be positive")
    heap: list[tuple[int, str]] = []
    records = 0
    bases = 0
    for identifier, lines in _records(path):
        records += 1
        bases += sum(len(line.strip()) for line in lines[1:])
        entry = (-_priority(identifier, seed), identifier)
        if len(heap) < count:
            heapq.heappush(heap, entry)
        elif entry > heap[0]:
            heapq.heapreplace(heap, entry)
    if records < count:
        raise ValueError(f"requested {count} records, but {path} contains only {records}")
    selected = {identifier for _, identifier in heap}
    if len(selected) != count:
        raise ValueError(f"{path}: duplicate FASTA identifiers prevent deterministic selection")
    return selected, SourceSummary(records=records, bases=bases)


def create_fixtures(
    source: Path,
    output_directory: Path,
    counts: tuple[int, ...],
    *,
    seed: str,
    force: bool = False,
) -> dict[str, object]:
    """Create nested FASTA/manifest fixtures and write non-versioned run metadata."""
    if not source.is_file():
        raise ValueError(f"source FASTA does not exist: {source}")
    if len(set(counts)) != len(counts):
        raise ValueError("record counts must be distinct")
    requested = tuple(sorted(counts))
    output_directory.mkdir(parents=True, exist_ok=True)
    destinations = {count: output_directory / f"real-contigs-{count}.fasta" for count in requested}
    manifests = {
        count: output_directory / f"real-contigs-{count}.samples.tsv" for count in requested
    }
    metadata_path = output_directory / "real-contigs-fixtures.json"
    output_paths = (*destinations.values(), *manifests.values(), metadata_path)
    existing = [path for path in output_paths if path.exists()]
    if existing and not force:
        rendered = ", ".join(str(path) for path in existing)
        raise ValueError(f"refusing to overwrite existing fixture(s): {rendered}")

    selected_largest, summary = select_identifiers(source, requested[-1], seed)
    ranked = sorted((_priority(identifier, seed), identifier) for identifier in selected_largest)
    selections = {count: {identifier for _, identifier in ranked[:count]} for count in requested}
    temporary_destinations = {
        count: destinations[count].with_suffix(".fasta.partial") for count in requested
    }
    written_records = {count: 0 for count in requested}
    written_bases = {count: 0 for count in requested}
    written_identifiers = {count: set() for count in requested}
    started = time.monotonic()
    try:
        with ExitStack() as stack:
            outputs = {
                count: stack.enter_context(
                    temporary_destinations[count].open("w", encoding="ascii")
                )
                for count in requested
            }
            for identifier, lines in _records(source):
                sequence_bases = sum(len(line.strip()) for line in lines[1:])
                raw_record = "".join(lines)
                for count in requested:
                    if identifier in selections[count]:
                        if identifier in written_identifiers[count]:
                            message = "duplicate FASTA identifiers prevent deterministic selection"
                            raise ValueError(f"{source}: {message}")
                        outputs[count].write(raw_record)
                        written_identifiers[count].add(identifier)
                        written_records[count] += 1
                        written_bases[count] += sequence_bases
        for count in requested:
            if written_records[count] != count:
                raise ValueError(f"{source}: selected {count} identifiers but wrote wrong count")
            os.replace(temporary_destinations[count], destinations[count])
    except BaseException:
        for path in temporary_destinations.values():
            path.unlink(missing_ok=True)
        raise
    elapsed_seconds = time.monotonic() - started
    for count in requested:
        sample = f"real_contigs_{count}"
        manifests[count].write_text(
            f"sample\tcontigs\n{sample}\t{destinations[count].name}\n", encoding="utf-8"
        )
    metadata = {
        "schema_version": 1,
        "source": {
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "records": summary.records,
            "bases": summary.bases,
        },
        "selection": {
            "method": "lowest blake2b identifier ranks; nested hash samples",
            "seed": seed,
            "records_requested": list(requested),
        },
        "fixtures": {
            str(count): {
                "fasta": destinations[count].name,
                "manifest": manifests[count].name,
                "records": written_records[count],
                "bases": written_bases[count],
            }
            for count in requested
        },
        "runtime": {
            "fixture_write_seconds": elapsed_seconds,
            "python": sys.version,
            "platform": platform.platform(),
        },
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def main() -> int:
    """Create fixtures from command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="plain or gzip FASTA source")
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--records", type=int, nargs="+", default=[10_000, 100_000, 1_000_000])
    parser.add_argument("--seed", default="contigger-real-fasta-scale-v1")
    parser.add_argument("--force", action="store_true", help="replace existing fixture paths")
    arguments = parser.parse_args()
    try:
        metadata = create_fixtures(
            arguments.source,
            arguments.output_directory,
            tuple(arguments.records),
            seed=arguments.seed,
            force=arguments.force,
        )
    except (OSError, UnicodeError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
