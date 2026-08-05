"""Stable exact sequence catalogue and provenance tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from contigger.catalogue import build_catalogue, catalogue_provenance, load_source_sequences
from contigger.exceptions import InputValidationError
from contigger.manifest import parse_manifest
from contigger.models import Orientation, SequenceRecord

DATASET = Path(__file__).parents[1] / "test_data"


def record(identifier: str, sequence: str, sample: str = "S") -> SequenceRecord:
    """Build a valid source record for catalogue tests."""
    return SequenceRecord(
        identifier=f"{sample}:{identifier}",
        source_sample=sample,
        original_identifier=identifier,
        description="",
        sequence=sequence,
        length=len(sequence),
    )


def test_exact_and_reverse_complement_deduplication_is_deterministic() -> None:
    records = (
        record("forward_b", "AACCGT", "B"),
        record("unique", "AAAAAC", "A"),
        record("reverse", "ACGGTT", "C"),
        record("forward_a", "AACCGT", "A"),
    )
    first = build_catalogue(records)
    assert first == build_catalogue(reversed(records))
    assert len(first.sequences) == 2
    duplicate_members = [
        member
        for member in first.members
        if member.original_identifier in {"forward_a", "forward_b", "reverse"}
    ]
    assert {member.orientation for member in duplicate_members} == {
        Orientation.FORWARD,
        Orientation.REVERSE,
    }
    assert sum(member.representative for member in duplicate_members) == 1


def test_palindrome_has_explicit_forward_orientation() -> None:
    catalogue = build_catalogue([record("palindrome", "ACGT")])
    assert catalogue.members[0].orientation is Orientation.FORWARD


def test_duplicate_source_identifier_is_rejected() -> None:
    duplicate = record("same", "ACGT")
    with pytest.raises(InputValidationError, match="duplicate source identifier"):
        build_catalogue([duplicate, duplicate])


def test_catalogue_provenance_retains_every_source() -> None:
    catalogue = build_catalogue([record("a", "AACCGT", "A"), record("b", "ACGGTT", "B")])
    provenance = catalogue_provenance(catalogue)
    assert len(provenance) == 2
    assert {item.source_contig for item in provenance} == {"a", "b"}
    assert {item.orientation for item in provenance} == {
        Orientation.FORWARD,
        Orientation.REVERSE,
    }
    assert all(item.relationship == "EXACT_MATCH" for item in provenance)


def test_pseudomonas_exact_truth_collapses_six_pairs() -> None:
    validation = parse_manifest(DATASET / "manifest.tsv")
    catalogue = build_catalogue(load_source_sequences(validation.samples))
    assert len(catalogue.members) == 90
    assert len(catalogue.sequences) == 84
    by_source = {member.original_identifier: member for member in catalogue.members}
    for prefix in ("exact", "reverse_exact"):
        for replicate in range(1, 4):
            left = by_source[f"{prefix}_r0{replicate}_a"]
            right = by_source[f"{prefix}_r0{replicate}_b"]
            assert left.catalogue_id == right.catalogue_id
            if prefix == "reverse_exact":
                assert left.orientation is not right.orientation
