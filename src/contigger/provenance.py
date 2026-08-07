"""Deterministic provenance records and TSV writer."""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from contigger.models import Orientation

PROVENANCE_COLUMNS = (
    "output_sequence",
    "source_sample",
    "source_contig",
    "relationship",
    "orientation",
    "source_start",
    "source_end",
    "output_start",
    "output_end",
    "identity",
    "disposition",
    "decision_reason",
    "evidence_mode",
)


@dataclass(frozen=True, slots=True)
class ProvenanceRecord:
    """Mapping from a source interval to an output interval."""

    output_sequence: str
    source_sample: str
    source_contig: str
    relationship: str
    orientation: Orientation
    source_start: int
    source_end: int
    output_start: int
    output_end: int
    identity: float | None
    disposition: str
    decision_reason: str
    evidence_mode: str = "none"


def write_provenance(path: Path, records: Iterable[ProvenanceRecord]) -> None:
    """Write provenance with fixed columns and deterministic row ordering."""
    ordered = sorted(
        records,
        key=lambda record: (
            record.output_sequence,
            record.output_start,
            record.source_sample,
            record.source_contig,
        ),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(PROVENANCE_COLUMNS)
        for record in ordered:
            writer.writerow(
                (
                    record.output_sequence,
                    record.source_sample,
                    record.source_contig,
                    record.relationship,
                    record.orientation.value,
                    record.source_start,
                    record.source_end,
                    record.output_start,
                    record.output_end,
                    "" if record.identity is None else f"{record.identity:.6f}",
                    record.disposition,
                    record.decision_reason,
                    record.evidence_mode,
                )
            )
