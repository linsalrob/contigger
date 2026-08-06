"""Selective alignment planning tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from contigger.alignment_planning import (
    execute_indexed_selective_alignments,
    execute_selective_alignments,
    plan_selective_alignments,
)
from contigger.exceptions import InputValidationError
from contigger.models import (
    AlignmentHit,
    CandidatePair,
    CatalogueSequence,
    Orientation,
    SequenceRecord,
)


class FakeAligner:
    """Minimal backend proving requests execute one pair at a time."""

    tool_name = "fake"
    tool_version = "1"
    last_command: tuple[str, ...] | None = None

    def build_index(self, targets: tuple[SequenceRecord, ...], index_path: object) -> object:
        return index_path

    def align(
        self, queries: tuple[SequenceRecord, ...], targets: tuple[SequenceRecord, ...]
    ) -> tuple[AlignmentHit, ...]:
        assert len(queries) == len(targets) == 1
        query, target = queries[0], targets[0]
        return (
            AlignmentHit(
                query.identifier,
                target.identifier,
                query.length,
                target.length,
                0,
                query.length,
                0,
                target.length,
                Orientation.FORWARD,
                query.length,
                query.length,
            ),
        )


class FakeIndexedAligner(FakeAligner):
    """Backend recording safe one-target query batches."""

    def __init__(self) -> None:
        self.batches: list[tuple[tuple[str, ...], tuple[str, ...]]] = []

    def align_indexed(
        self,
        queries: tuple[SequenceRecord, ...],
        targets: tuple[SequenceRecord, ...],
        index_path: object,
    ) -> tuple[AlignmentHit, ...]:
        self.batches.append(
            (
                tuple(record.identifier for record in queries),
                tuple(record.identifier for record in targets),
            )
        )
        return tuple(
            AlignmentHit(
                query.identifier,
                targets[0].identifier,
                query.length,
                targets[0].length,
                0,
                query.length,
                0,
                targets[0].length,
                Orientation.FORWARD,
                query.length,
                query.length,
            )
            for query in queries
        )


def sequence(identifier: str) -> CatalogueSequence:
    """Build a small catalogue sequence."""
    return CatalogueSequence(identifier, "ACGT", 4, identifier, identifier)


def test_requests_resolve_exact_candidate_pairs_deterministically() -> None:
    candidate = CandidatePair("a", "b", 7)
    requests = plan_selective_alignments([sequence("b"), sequence("a")], [candidate])
    assert len(requests) == 1
    assert requests[0].query.identifier == "a"
    assert requests[0].target.identifier == "b"
    assert requests[0].candidate is candidate


def test_unknown_sequence_and_noncanonical_pair_are_rejected() -> None:
    with pytest.raises(InputValidationError, match="unknown catalogue"):
        plan_selective_alignments([sequence("a")], [CandidatePair("a", "b", 2)])
    with pytest.raises(InputValidationError, match="canonical identifier order"):
        plan_selective_alignments([sequence("a"), sequence("b")], [CandidatePair("b", "a", 2)])


def test_executor_calls_backend_for_only_the_planned_pair() -> None:
    requests = plan_selective_alignments(
        [sequence("a"), sequence("b")], [CandidatePair("a", "b", 2)]
    )
    hits = execute_selective_alignments(requests, FakeAligner())
    assert [(hit.query_id, hit.target_id) for hit in hits] == [("a", "b")]


def test_indexed_executor_batches_queries_only_for_their_approved_target(tmp_path: Path) -> None:
    requests = plan_selective_alignments(
        [sequence("a"), sequence("b"), sequence("c")],
        [CandidatePair("a", "c", 2), CandidatePair("b", "c", 2), CandidatePair("a", "b", 2)],
    )
    aligner = FakeIndexedAligner()
    hits = execute_indexed_selective_alignments(requests, aligner, tmp_path)
    assert aligner.batches == [(("a",), ("b",)), (("a", "b"), ("c",))]
    assert [(hit.query_id, hit.target_id) for hit in hits] == [
        ("a", "b"),
        ("a", "c"),
        ("b", "c"),
    ]
