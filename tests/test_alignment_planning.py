"""Selective alignment planning tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from contigger.alignment_planning import (
    execute_indexed_selective_alignments,
    execute_selective_alignments,
    iter_indexed_selective_alignment_batches,
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


class OffBatchAligner(FakeIndexedAligner):
    """Backend returning another planned pair from the wrong target batch."""

    def align_indexed(
        self,
        queries: tuple[SequenceRecord, ...],
        targets: tuple[SequenceRecord, ...],
        index_path: object,
    ) -> tuple[AlignmentHit, ...]:
        if targets[0].identifier == "b":
            return (
                AlignmentHit(
                    "a",
                    "c",
                    4,
                    4,
                    0,
                    4,
                    0,
                    4,
                    Orientation.FORWARD,
                    4,
                    4,
                ),
            )
        return super().align_indexed(queries, targets, index_path)


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


def test_indexed_executor_reuses_one_target_index_across_query_batches(tmp_path: Path) -> None:
    """Large target groups are split without introducing cross-target alignments."""
    requests = plan_selective_alignments(
        [sequence("a"), sequence("b"), sequence("c")],
        [CandidatePair("a", "c", 2), CandidatePair("b", "c", 2)],
    )
    aligner = FakeIndexedAligner()
    hits = execute_indexed_selective_alignments(
        requests, aligner, tmp_path, max_queries_per_batch=1
    )
    assert aligner.batches == [(("a",), ("c",)), (("b",), ("c",))]
    assert [(hit.query_id, hit.target_id) for hit in hits] == [("a", "c"), ("b", "c")]


def test_indexed_batch_iterator_does_not_accumulate_other_batches(tmp_path: Path) -> None:
    """The streaming API yields one deterministic target/query batch at a time."""
    requests = plan_selective_alignments(
        [sequence("a"), sequence("b"), sequence("c")],
        [CandidatePair("a", "c", 2), CandidatePair("b", "c", 2)],
    )
    aligner = FakeIndexedAligner()
    batches = iter_indexed_selective_alignment_batches(
        requests, aligner, tmp_path, max_queries_per_batch=1
    )

    first = next(batches)
    assert [(hit.query_id, hit.target_id) for hit in first] == [("a", "c")]
    assert aligner.batches == [(("a",), ("c",))]
    assert [(hit.query_id, hit.target_id) for hit in next(batches)] == [("b", "c")]
    with pytest.raises(StopIteration):
        next(batches)


def test_indexed_batch_iterator_reports_bounded_progress(tmp_path: Path) -> None:
    """Progress callbacks receive completed batches and observations."""
    requests = plan_selective_alignments(
        [sequence("a"), sequence("b"), sequence("c")],
        [CandidatePair("a", "c", 2), CandidatePair("b", "c", 2)],
    )
    progress: list[tuple[int, int, int]] = []
    list(
        iter_indexed_selective_alignment_batches(
            requests,
            FakeIndexedAligner(),
            tmp_path,
            max_queries_per_batch=1,
            progress_callback=lambda completed, total, observations: progress.append(
                (completed, total, observations)
            ),
        )
    )

    assert progress == [(1, 2, 1), (2, 2, 1)]


def test_indexed_executor_rejects_pair_from_another_planned_batch(tmp_path: Path) -> None:
    requests = plan_selective_alignments(
        [sequence("a"), sequence("b"), sequence("c")],
        [CandidatePair("a", "b", 2), CandidatePair("a", "c", 2)],
    )
    with pytest.raises(InputValidationError, match="current selective batch"):
        execute_indexed_selective_alignments(requests, OffBatchAligner(), tmp_path)
