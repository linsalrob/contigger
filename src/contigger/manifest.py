"""Line-aware TSV sample manifest parsing and validation."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from contigger.exceptions import ManifestError
from contigger.fasta import validate_fasta
from contigger.models import SampleInput

REQUIRED_COLUMNS = frozenset({"sample", "contigs"})
OPTIONAL_COLUMNS = frozenset({"bam", "technology", "assembly_graph"})


@dataclass(frozen=True, slots=True)
class ManifestValidation:
    """Validated, deterministically ordered samples and non-fatal warnings."""

    samples: tuple[SampleInput, ...]
    warnings: tuple[str, ...] = ()


def parse_manifest(path: Path, *, check_files: bool = True) -> ManifestValidation:
    """Parse a TSV manifest, resolving file paths relative to its directory."""
    try:
        handle = path.open(encoding="utf-8", newline="")
    except OSError as error:
        raise ManifestError(f"cannot read manifest {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ManifestError(f"{path}:1: manifest has no header")
        columns = [column.strip() for column in reader.fieldnames]
        if len(columns) != len(set(columns)):
            raise ManifestError(f"{path}:1: duplicate manifest column")
        missing = sorted(REQUIRED_COLUMNS - set(columns))
        if missing:
            raise ManifestError(f"{path}:1: missing required column(s): {', '.join(missing)}")
        reader.fieldnames = columns
        samples: list[SampleInput] = []
        seen: set[str] = set()
        warnings: list[str] = []
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise ManifestError(f"{path}:{line_number}: too many tab-separated fields")
            row = {key: (value or "").strip() for key, value in raw_row.items()}
            sample = row["sample"]
            contigs_value = row["contigs"]
            if not sample or not contigs_value:
                raise ManifestError(f"{path}:{line_number}: sample and contigs values are required")
            if sample in seen:
                raise ManifestError(f"{path}:{line_number}: duplicate sample identifier {sample!r}")
            seen.add(sample)
            contigs = _resolve(path.parent, contigs_value)
            bam = _optional_path(path.parent, row.get("bam", ""))
            graph = _optional_path(path.parent, row.get("assembly_graph", ""))
            metadata = {
                key: value
                for key, value in row.items()
                if key not in REQUIRED_COLUMNS | OPTIONAL_COLUMNS
            }
            item = SampleInput(
                sample=sample,
                contigs=contigs,
                bam=bam,
                technology=row.get("technology") or None,
                assembly_graph=graph,
                metadata=metadata,
            )
            if check_files:
                _validate_sample_files(item, path, line_number, warnings)
            samples.append(item)
    if not samples:
        raise ManifestError(f"{path}: manifest contains no sample rows")
    return ManifestValidation(tuple(sorted(samples, key=lambda item: item.sample)), tuple(warnings))


def _resolve(base: Path, value: str) -> Path:
    candidate = Path(value).expanduser()
    return (base / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()


def _optional_path(base: Path, value: str) -> Path | None:
    return _resolve(base, value) if value else None


def _validate_sample_files(
    sample: SampleInput, manifest: Path, line_number: int, warnings: list[str]
) -> None:
    if not sample.contigs.is_file():
        raise ManifestError(
            f"{manifest}:{line_number}: contig FASTA does not exist: {sample.contigs}"
        )
    validate_fasta(sample.contigs, sample.sample)
    optional_inputs = (("BAM/CRAM", sample.bam), ("assembly graph", sample.assembly_graph))
    for label, optional_path in optional_inputs:
        if optional_path is not None and not optional_path.is_file():
            raise ManifestError(
                f"{manifest}:{line_number}: {label} file does not exist: {optional_path}"
            )
    if sample.bam is not None and not _has_alignment_index(sample.bam):
        warnings.append(f"sample {sample.sample}: BAM/CRAM lacks an adjacent index: {sample.bam}")


def _has_alignment_index(path: Path) -> bool:
    candidates = [Path(f"{path}.bai"), Path(f"{path}.crai")]
    if path.suffix.lower() == ".bam":
        candidates.append(path.with_suffix(".bai"))
    elif path.suffix.lower() == ".cram":
        candidates.append(path.with_suffix(".crai"))
    return any(candidate.is_file() for candidate in candidates)
