"""Canonical positional-minimiser candidate tests."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

import contigger.minimisers as minimisers
from contigger.benchmark import MERGE_LIKE, load_truth
from contigger.catalogue import build_catalogue, load_source_sequences
from contigger.exceptions import ConfigurationError
from contigger.manifest import parse_manifest
from contigger.minimisers import (
    generate_candidates,
    generate_candidates_with_metrics,
    sequence_minimisers,
)
from contigger.models import CatalogueSequence, Orientation
from contigger.utilities.sequences import reverse_complement

DATASET = Path(__file__).parents[1] / "test_data"


def sequence(identifier: str, bases: str) -> CatalogueSequence:
    """Build a catalogue sequence with an irrelevant test digest."""
    return CatalogueSequence(identifier, bases, len(bases), identifier, identifier)


def reference_minimisers(
    item: CatalogueSequence, *, kmer_size: int, window_size: int
) -> tuple[tuple[int, int, Orientation, str], ...]:
    """Implement the straightforward sliding-window definition for regression tests."""
    kmers: list[tuple[int, int, Orientation, str] | None] = []
    for position in range(item.length - kmer_size + 1):
        forward = item.sequence[position : position + kmer_size]
        if set(forward) - set("ACGT"):
            kmers.append(None)
            continue
        reverse = reverse_complement(forward)
        canonical = min(forward, reverse)
        orientation = Orientation.FORWARD if forward == canonical else Orientation.REVERSE
        value = int.from_bytes(hashlib.sha256(canonical.encode("ascii")).digest()[:8], "big")
        kmers.append((value, position, orientation, canonical))
    selected: set[tuple[int, int, Orientation, str]] = set()
    effective_window = min(window_size, len(kmers))
    for start in range(len(kmers) - effective_window + 1):
        valid = [value for value in kmers[start : start + effective_window] if value is not None]
        if valid:
            minimum = min(value[0] for value in valid)
            selected.update(value for value in valid if value[0] == minimum)
    return tuple(sorted(selected, key=lambda value: (value[1], value[0], value[2].value, value[3])))


def test_minimisers_are_deterministic_and_strand_explicit() -> None:
    item = sequence("a", "AACCGTTA")
    first = sequence_minimisers(item, kmer_size=3, window_size=2)
    assert first == sequence_minimisers(item, kmer_size=3, window_size=2)
    assert first == tuple(sorted(first, key=lambda value: value.position))
    assert {observation.orientation for observation in first} <= {
        Orientation.FORWARD,
        Orientation.REVERSE,
    }


def test_monotonic_minimiser_selection_matches_window_definition() -> None:
    for bases in ("AACCGTTA", "AAAAAAA", "ACGTNNNACGT", "GATTACAGATTACA"):
        item = sequence("a", bases)
        observed = tuple(
            (value.value, value.position, value.orientation, value.kmer)
            for value in sequence_minimisers(item, kmer_size=3, window_size=3)
        )
        assert observed == reference_minimisers(item, kmer_size=3, window_size=3)


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


@pytest.mark.parametrize("limit", [0, -1])
def test_seed_pair_limit_must_be_positive(limit: int) -> None:
    with pytest.raises(ConfigurationError, match="seed-pair"):
        generate_candidates(
            [sequence("a", "AACCGGTT")],
            kmer_size=3,
            window_size=2,
            min_shared_minimisers=1,
            max_minimiser_frequency=20,
            terminal_band=2,
            max_seed_pair_observations=limit,
        )


def test_candidate_shard_count_is_bounded() -> None:
    with pytest.raises(ConfigurationError, match="shard count"):
        generate_candidates(
            [sequence("a", "AACCGGTT")],
            kmer_size=3,
            window_size=2,
            min_shared_minimisers=1,
            max_minimiser_frequency=20,
            terminal_band=2,
            candidate_shards=65,
        )


def test_candidate_results_are_identical_across_shard_counts() -> None:
    shared = "ACGTCAGTACGATCGTACGA"
    sequences = [
        sequence("a", "TTGCAACGTT" + shared),
        sequence("b", shared + "GGCATTAACC"),
        sequence("c", shared + "TTAGGCCAAT"),
    ]
    options = {
        "kmer_size": 5,
        "window_size": 3,
        "min_shared_minimisers": 2,
        "max_minimiser_frequency": 20,
        "terminal_band": 8,
    }
    assert generate_candidates(sequences, candidate_shards=1, **options) == generate_candidates(
        sequences, candidate_shards=3, **options
    )


def test_external_seed_sort_chunks_preserve_candidate_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Bounded external sorting must not change positional candidate evidence."""
    shared = "ACGTCAGTACGATCGTACGA"
    sequences = [
        sequence("a", "TTGCAACGTT" + shared),
        sequence("b", shared + "GGCATTAACC"),
        sequence("c", shared + "TTAGGCCAAT"),
    ]
    options = {
        "kmer_size": 5,
        "window_size": 3,
        "min_shared_minimisers": 2,
        "max_minimiser_frequency": 20,
        "terminal_band": 8,
        "candidate_shards": 1,
    }
    expected = generate_candidates(sequences, **options)
    monkeypatch.setattr(minimisers, "_SEED_SORT_CHUNK_LINES", 1)
    monkeypatch.setattr(minimisers, "_SEED_SORT_FAN_IN", 2)
    assert generate_candidates(sequences, **options) == expected


def test_seed_sort_fan_in_reserves_pair_writer_descriptors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """External sorting leaves headroom for its bounded pair-output writer pool."""
    monkeypatch.setattr(minimisers.resource, "getrlimit", lambda _resource: (64, 64))
    assert minimisers._seed_sort_fan_in() == 39


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
