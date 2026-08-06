"""Validated selective-alignment requests derived from candidate evidence."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

from contigger.aligners.base import Aligner, IndexedAligner
from contigger.exceptions import InputValidationError
from contigger.models import (
    AlignmentHit,
    AlignmentRequest,
    CandidatePair,
    CatalogueSequence,
    SequenceRecord,
)
from contigger.relationships import alignment_sort_key


def plan_selective_alignments(
    sequences: Iterable[CatalogueSequence], candidates: Iterable[CandidatePair]
) -> tuple[AlignmentRequest, ...]:
    """Resolve candidate identifiers to exact per-pair alignment requests."""
    sequence_items = tuple(sequences)
    lookup = {item.identifier: item for item in sequence_items}
    if len(lookup) != len(sequence_items):
        raise InputValidationError("selective alignment requires unique catalogue identifiers")
    requests: list[AlignmentRequest] = []
    seen: set[tuple[str, str]] = set()
    for candidate in sorted(candidates, key=lambda item: (item.query_id, item.target_id)):
        pair = (candidate.query_id, candidate.target_id)
        if pair[0] >= pair[1]:
            raise InputValidationError("alignment candidates must use canonical identifier order")
        if pair in seen:
            raise InputValidationError(f"duplicate selective-alignment candidate: {pair}")
        seen.add(pair)
        try:
            query = lookup[pair[0]]
            target = lookup[pair[1]]
        except KeyError as error:
            raise InputValidationError(
                f"candidate references unknown catalogue sequence: {error.args[0]}"
            ) from error
        requests.append(
            AlignmentRequest(
                query=_as_record(query),
                target=_as_record(target),
                candidate=candidate,
            )
        )
    return tuple(requests)


def execute_selective_alignments(
    requests: Iterable[AlignmentRequest], aligner: Aligner
) -> tuple[AlignmentHit, ...]:
    """Execute only planned pairs through a replaceable alignment backend.

    Each backend call receives one query and one target, preventing accidental
    all-v-all expansion. Returned identifiers must match that exact ordered request.
    """
    hits: list[AlignmentHit] = []
    for request in sorted(
        requests, key=lambda item: (item.query.identifier, item.target.identifier)
    ):
        for hit in aligner.align((request.query,), (request.target,)):
            expected = (request.query.identifier, request.target.identifier)
            observed = (hit.query_id, hit.target_id)
            if observed != expected:
                raise InputValidationError(
                    "alignment backend returned identifiers outside selective request: "
                    f"expected {expected}, observed {observed}"
                )
            hits.append(hit)
    return tuple(sorted(hits, key=alignment_sort_key))


def execute_indexed_selective_alignments(
    requests: Iterable[AlignmentRequest],
    aligner: IndexedAligner,
    index_directory: Path,
) -> tuple[AlignmentHit, ...]:
    """Batch approved queries by one target and reuse validated target indexes.

    A batch never contains multiple targets, so minimap2 cannot introduce
    unrequested query-target combinations through an all-v-all expansion.
    """
    request_items = tuple(
        sorted(requests, key=lambda item: (item.target.identifier, item.query.identifier))
    )
    expected: set[tuple[str, str]] = set()
    grouped: dict[str, list[AlignmentRequest]] = defaultdict(list)
    for request in request_items:
        pair = (request.query.identifier, request.target.identifier)
        if pair in expected:
            raise InputValidationError(f"duplicate selective-alignment request: {pair}")
        expected.add(pair)
        grouped[request.target.identifier].append(request)

    hits: list[AlignmentHit] = []
    for target_id in sorted(grouped):
        group = grouped[target_id]
        target = group[0].target
        if any(request.target != target for request in group):
            raise InputValidationError(
                f"target identifier collision in selective alignment batch: {target_id}"
            )
        index_name = hashlib.sha256(target_id.encode("utf-8")).hexdigest()
        index_path = index_directory / f"target-{index_name}.mmi"
        aligner.build_index((target,), index_path)
        queries = tuple(request.query for request in group)
        for hit in aligner.align_indexed(queries, (target,), index_path):
            observed = (hit.query_id, hit.target_id)
            if observed not in expected:
                raise InputValidationError(
                    "alignment backend returned identifiers outside selective requests: "
                    f"observed {observed}"
                )
            hits.append(hit)
    return tuple(sorted(hits, key=alignment_sort_key))


def _as_record(sequence: CatalogueSequence) -> SequenceRecord:
    return SequenceRecord(
        identifier=sequence.identifier,
        source_sample="",
        original_identifier=sequence.identifier,
        description="canonical catalogue sequence",
        sequence=sequence.sequence,
        length=sequence.length,
    )
