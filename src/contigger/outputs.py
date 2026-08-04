"""Centralised output naming for implemented and planned artifacts."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class OutputPaths:
    """All names derived consistently from one output prefix."""

    fasta: Path
    provenance: Path
    relationships: Path
    ambiguous: Path
    gfa: Path
    stats: Path
    variants: Path
    join_support: Path
    consensus_vcf: Path
    low_confidence_bed: Path


def output_paths(prefix: Path) -> OutputPaths:
    """Resolve all planned artifact names without creating files."""
    stem = str(prefix)
    return OutputPaths(
        fasta=Path(f"{stem}.fasta"),
        provenance=Path(f"{stem}.provenance.tsv"),
        relationships=Path(f"{stem}.relationships.tsv"),
        ambiguous=Path(f"{stem}.ambiguous.tsv"),
        gfa=Path(f"{stem}.gfa"),
        stats=Path(f"{stem}.stats.json"),
        variants=Path(f"{stem}.variants.tsv"),
        join_support=Path(f"{stem}.join_support.tsv"),
        consensus_vcf=Path(f"{stem}.consensus.vcf"),
        low_confidence_bed=Path(f"{stem}.low_confidence.bed"),
    )
