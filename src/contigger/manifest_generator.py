"""Discover assembly sidecars and write a valid Contigger manifest."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

FASTA_SUFFIXES = (
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fas.gz",
    ".fsa.gz",
    ".fasta",
    ".fa",
    ".fna",
    ".fas",
    ".fsa",
)
GRAPH_SUFFIXES = (".gfa.gz", ".gfa")
ALIGNMENT_SUFFIXES = (".bam", ".cram")


@dataclass(frozen=True, slots=True)
class ManifestRow:
    """One discovered assembly and its optional sidecars."""

    sample: str
    contigs: Path
    bam: Path | None
    assembly_graph: Path | None
    index: Path | None


def _stem(path: Path, suffixes: tuple[str, ...]) -> str | None:
    """Remove one recognized suffix from a filename."""
    lowered = path.name.casefold()
    for suffix in suffixes:
        if lowered.endswith(suffix):
            return path.name[: -len(suffix)]
    return None


def sample_name(path: Path) -> str | None:
    """Return a likely sample name for a FASTA filename."""
    return _stem(path, FASTA_SUFFIXES)


def _index(paths: list[Path], suffixes: tuple[str, ...]) -> dict[str, list[Path]]:
    result: dict[str, list[Path]] = {}
    for path in paths:
        stem = _stem(path, suffixes)
        if stem is not None:
            result.setdefault(stem.casefold(), []).append(path)
    return result


def _match(index: dict[str, list[Path]], stem: str, label: str) -> Path | None:
    matches = index.get(stem.casefold(), [])
    if len(matches) > 1:
        names = ", ".join(str(item) for item in matches)
        raise ValueError(f"multiple {label} files match sample {stem!r}: {names}")
    return matches[0] if matches else None


def _alignment_index(path: Path) -> Path | None:
    if path.suffix.casefold() == ".bam":
        candidates = (Path(f"{path}.bai"), path.with_suffix(".bai"))
    else:
        candidates = (Path(f"{path}.crai"), path.with_suffix(".crai"))
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def discover(directory: Path, *, recursive: bool = True) -> tuple[ManifestRow, ...]:
    """Find FASTA files and matching GFA/BAM/CRAM sidecars."""
    if not directory.is_dir():
        raise ValueError(f"input directory does not exist or is not a directory: {directory}")
    iterator = directory.rglob("*") if recursive else directory.glob("*")
    files = sorted((item for item in iterator if item.is_file()), key=str)
    assemblies = [
        (item, stem) for item in files if (stem := sample_name(item)) is not None and stem
    ]
    if not assemblies:
        raise ValueError(f"no FASTA files found in {directory}")
    graphs = _index(files, GRAPH_SUFFIXES)
    alignments = _index(files, ALIGNMENT_SUFFIXES)
    rows: list[ManifestRow] = []
    seen: set[str] = set()
    for contigs, stem in sorted(assemblies, key=lambda item: (item[1].casefold(), str(item[0]))):
        if stem.casefold() in seen:
            raise ValueError(f"multiple FASTA files produce sample name {stem!r}")
        seen.add(stem.casefold())
        bam = _match(alignments, stem, "alignment")
        rows.append(
            ManifestRow(
                stem,
                contigs,
                bam,
                _match(graphs, stem, "assembly graph"),
                _alignment_index(bam) if bam else None,
            )
        )
    return tuple(rows)


def write_manifest(rows: tuple[ManifestRow, ...], output: Path) -> tuple[str, ...]:
    """Write a tab-separated manifest and return non-fatal warnings."""
    output.parent.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    def relative(path: Path) -> str:
        return (
            path.relative_to(output.parent).as_posix()
            if path.is_relative_to(output.parent)
            else path.as_posix()
        )

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "contigs", "bam", "technology", "assembly_graph"))
        for row in rows:
            if row.bam and row.index is None:
                warnings.append(
                    f"sample {row.sample}: no adjacent BAM/CRAM index found for {row.bam}"
                )
            writer.writerow(
                (
                    row.sample,
                    relative(row.contigs),
                    relative(row.bam) if row.bam else "",
                    "",
                    relative(row.assembly_graph) if row.assembly_graph else "",
                )
            )
    return tuple(warnings)


def create_manifest(
    directory: Path, output: Path, *, recursive: bool = True
) -> tuple[tuple[ManifestRow, ...], tuple[str, ...]]:
    """Discover inputs and write a manifest."""
    rows = discover(directory.resolve(), recursive=recursive)
    return rows, write_manifest(rows, output.resolve())
