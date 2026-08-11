"""Canonical positional-minimiser candidate tests."""

from __future__ import annotations

from pathlib import Path

from contigger.benchmark import MERGE_LIKE, load_truth
from contigger.catalogue import build_catalogue, load_source_sequences
from contigger.manifest import parse_manifest
from contigger.minimisers import (
    generate_candidates,
    generate_candidates_with_metrics,
    sequence_minimisers,
)
from contigger.models import CatalogueSequence, Orientation

DATASET = Path(__file__).parents[1] / "test_data"


def sequence(identifier: str, bases: str) -> CatalogueSequence:
    """Build a catalogue sequence with an irrelevant test digest."""
    return CatalogueSequence(identifier, bases, len(bases), identifier, identifier)


def test_minimisers_are_deterministic_and_strand_explicit() -> None:
    item = sequence("a", "AACCGTTA")
    first = sequence_minimisers(item, kmer_size=3, window_size=2)
    assert first == sequence_minimisers(item, kmer_size=3, window_size=2)
    assert first == tuple(sorted(first, key=lambda value: value.position))
    assert {observation.orientation for observation in first} <= {
        Orientation.FORWARD,
        Orientation.REVERSE,
    }


def test_ambiguous_kmers_are_not_seed_evidence() -> None:
    observations = sequence_minimisers(sequence("a", "AAANNNTTT"), kmer_size=3, window_size=2)
    assert all("N" not in item.kmer for item in observations)


def test_terminal_overlap_emits_candidate_with_position_evidence() -> None:
    shared = "ACGTCAGTACGATCGTACGA"
    left = sequence("a", "TTGCAACGTT" + shared)
    right = sequence("b", shared + "GGCATTAACC")
    candidates = generate_candidates(
        [right, left],
        kmer_size=5,
        window_size=3,
        min_shared_minimisers=2,
        max_minimiser_frequency=20,
        terminal_band=8,
    )
    assert len(candidates) == 1
    candidate = candidates[0]
    assert candidate.query_positions
    assert candidate.target_positions
    assert any("SUFFIX_TO_TARGET_PREFIX" in item for item in candidate.terminal_topologies)


def test_candidate_metrics_describe_retained_seed_pressure() -> None:
    shared = "ACGTCAGTACGATCGTACGA"
    candidates, metrics = generate_candidates_with_metrics(
        [sequence("a", "TTGCAACGTT" + shared), sequence("b", shared + "GGCATTAACC")],
        kmer_size=5,
        window_size=3,
        min_shared_minimisers=2,
        max_minimiser_frequency=20,
        terminal_band=8,
    )
    assert len(candidates) == metrics.candidate_pairs == 1
    assert metrics.input_sequences == 2
    assert metrics.input_bases == 60
    assert metrics.retained_observations <= metrics.minimiser_observations
    assert metrics.unique_minimisers > 0


def test_internal_similarity_does_not_emit_candidate() -> None:
    shared = "ACGTCAGTACGATCGTACGA"
    left = sequence("a", "A" * 15 + shared + "C" * 15)
    right = sequence("b", "G" * 15 + shared + "T" * 15)
    assert not generate_candidates(
        [left, right],
        kmer_size=5,
        window_size=3,
        min_shared_minimisers=2,
        max_minimiser_frequency=20,
        terminal_band=5,
    )


def test_frequent_minimisers_are_suppressed() -> None:
    sequences = [sequence(f"s{index}", "A" * 30) for index in range(3)]
    assert not generate_candidates(
        sequences,
        kmer_size=5,
        window_size=3,
        min_shared_minimisers=1,
        max_minimiser_frequency=2,
        terminal_band=5,
    )


def test_pseudomonas_valid_pairwise_cases_reach_selective_alignment() -> None:
    validation = parse_manifest(DATASET / "manifest.tsv")
    catalogue = build_catalogue(load_source_sequences(validation.samples))
    candidates = generate_candidates(
        catalogue.sequences,
        kmer_size=21,
        window_size=10,
        min_shared_minimisers=5,
        max_minimiser_frequency=100,
        terminal_band=1000,
    )
    candidate_pairs = {(item.query_id, item.target_id) for item in candidates}
    source_to_catalogue = {
        member.original_identifier: member.catalogue_id for member in catalogue.members
    }
    truth = load_truth(DATASET / "expected" / "expected_relationships.tsv")
    expected_cases = {
        item.case_id
        for item in truth
        if item.merge_allowed
        and item.relationship_type in MERGE_LIKE
        and item.relationship_type.value != "EXACT_MATCH"
    }
    observed_cases = {
        item.case_id
        for item in truth
        if tuple(
            sorted(
                (
                    source_to_catalogue[item.query_id],
                    source_to_catalogue[item.target_id],
                )
            )
        )
        in candidate_pairs
    }
    assert expected_cases <= observed_cases
    assert len(candidates) == 61
