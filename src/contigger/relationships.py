"""Conservative single-hit and complete ordered-pair classification."""

from collections.abc import Iterable, Iterator
from dataclasses import replace

from contigger.exceptions import InputValidationError
from contigger.models import (
    AlignmentHit,
    PairRelationship,
    RejectedAlignment,
    Relationship,
    RelationshipType,
    RunConfig,
)


def classify_alignment(hit: AlignmentHit, config: RunConfig) -> Relationship:
    """Classify one alignment using identity, coverage, and terminal geometry.

    This initial classifier intentionally does not resolve competing alignments or
    repeat signals. Those inputs must be added before production graph construction.
    """
    reasons: list[str] = []
    if hit.alignment_block_length == 0:
        return _relationship(hit, RelationshipType.NO_RELATIONSHIP, "rejected", ("zero length",))
    if hit.identity < config.identity:
        return _relationship(
            hit,
            RelationshipType.NO_RELATIONSHIP,
            "rejected",
            ("below identity threshold",),
        )

    tolerance = config.end_tolerance
    query_full = hit.query_start_distance <= tolerance and hit.query_end_distance <= tolerance
    target_full = hit.target_start_distance <= tolerance and hit.target_end_distance <= tolerance
    if (
        query_full
        and target_full
        and hit.query_coverage >= config.containment_coverage
        and hit.target_coverage >= config.containment_coverage
        and hit.identity == 1.0
        and hit.query_length == hit.target_length
    ):
        return _relationship(hit, RelationshipType.EXACT_MATCH, "confident", ("full exact span",))

    if (
        hit.alignment_block_length >= config.min_containment
        and query_full
        and hit.query_coverage >= config.containment_coverage
        and hit.query_length <= hit.target_length
    ):
        return _relationship(
            hit,
            RelationshipType.QUERY_CONTAINED_IN_TARGET,
            "candidate",
            ("query covered end-to-end",),
        )
    if (
        hit.alignment_block_length >= config.min_containment
        and target_full
        and hit.target_coverage >= config.containment_coverage
        and hit.target_length <= hit.query_length
    ):
        return _relationship(
            hit,
            RelationshipType.TARGET_CONTAINED_IN_QUERY,
            "candidate",
            ("target covered end-to-end",),
        )

    if hit.alignment_block_length < config.min_overlap:
        return _relationship(
            hit, RelationshipType.NO_RELATIONSHIP, "rejected", ("below minimum overlap",)
        )

    # PAF target coordinates always refer to the forward target. Swap their semantic
    # ends when the query aligns to the reverse-complement target orientation.
    oriented_target_start = (
        hit.target_start_distance if hit.orientation.value == "+" else hit.target_end_distance
    )
    oriented_target_end = (
        hit.target_end_distance if hit.orientation.value == "+" else hit.target_start_distance
    )
    query_suffix_target_prefix = (
        hit.query_end_distance <= tolerance and oriented_target_start <= tolerance
    )
    target_suffix_query_prefix = (
        hit.query_start_distance <= tolerance and oriented_target_end <= tolerance
    )
    if query_suffix_target_prefix and target_suffix_query_prefix:
        return _relationship(
            hit,
            RelationshipType.AMBIGUOUS_OVERLAP,
            "ambiguous",
            ("alignment satisfies competing terminal topologies",),
        )
    if query_suffix_target_prefix:
        return _relationship(
            hit,
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
            "candidate",
            ("compatible terminal geometry",),
        )
    if target_suffix_query_prefix:
        return _relationship(
            hit,
            RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX,
            "candidate",
            ("compatible terminal geometry",),
        )
    reasons.append("alignment lacks compatible terminal geometry")
    return _relationship(hit, RelationshipType.NO_RELATIONSHIP, "rejected", tuple(reasons))


def relationship_sort_key(relationship: Relationship) -> tuple[str, str, str, str]:
    """Return the documented deterministic output ordering key."""
    return (
        relationship.query_id,
        relationship.target_id,
        relationship.relationship_type.value,
        relationship.orientation.value,
    )


def alignment_sort_key(hit: AlignmentHit) -> tuple[object, ...]:
    """Return a total, deterministic ordering key for an alignment record."""
    return (
        hit.query_id,
        hit.target_id,
        hit.orientation.value,
        hit.query_length,
        hit.target_length,
        hit.query_start,
        hit.query_end,
        hit.target_start,
        hit.target_end,
        hit.matching_bases,
        hit.alignment_block_length,
        -1 if hit.mapping_quality is None else hit.mapping_quality,
        -1 if hit.alignment_score is None else hit.alignment_score,
        -1 if hit.supporting_seeds is None else hit.supporting_seeds,
        -1 if hit.chaining_score is None else hit.chaining_score,
        -1 if hit.secondary_chaining_score is None else hit.secondary_chaining_score,
        "" if hit.alignment_type is None else hit.alignment_type.value,
    )


def group_ordered_pairs(hits: Iterable[AlignmentHit]) -> Iterator[tuple[AlignmentHit, ...]]:
    """Group a complete hit stream by ordered query-target pair deterministically."""
    groups: dict[tuple[str, str], list[AlignmentHit]] = {}
    for hit in hits:
        groups.setdefault((hit.query_id, hit.target_id), []).append(hit)
    for pair in sorted(groups):
        yield tuple(sorted(groups[pair], key=alignment_sort_key))


def classify_pair(hits: Iterable[AlignmentHit], config: RunConfig) -> PairRelationship:
    """Classify all distinct alignments for one ordered query-target pair.

    Accepted hits may collapse only when relationship type, orientation, and every
    placement coordinate differ by no more than ``end_tolerance``. Any other
    competing accepted evidence is retained and reported as ambiguous; scores and
    primary/secondary labels never elect a winner.
    """
    ordered = tuple(sorted(set(hits), key=alignment_sort_key))
    if not ordered:
        raise InputValidationError("pair classification requires at least one alignment")
    pair = (ordered[0].query_id, ordered[0].target_id)
    if any((hit.query_id, hit.target_id) != pair for hit in ordered):
        raise InputValidationError("pair classification received multiple ordered pairs")

    accepted: list[tuple[AlignmentHit, Relationship]] = []
    rejected: list[RejectedAlignment] = []
    for hit in ordered:
        relationship = classify_alignment(hit, config)
        if relationship.relationship_type is RelationshipType.NO_RELATIONSHIP:
            rejected.append(RejectedAlignment(hit, relationship))
        else:
            accepted.append((hit, relationship))

    if not accepted:
        rejected_representative = ordered[0]
        decision = classify_alignment(rejected_representative, config)
        reasons = tuple(
            sorted({reason for item in rejected for reason in item.relationship.reasons})
        )
        return PairRelationship(
            replace(decision, reasons=reasons), rejected_representative, (), tuple(rejected)
        )

    representative_hit, representative_relationship = accepted[0]
    incompatibilities = _incompatibilities(accepted, config.end_tolerance)
    if incompatibilities:
        decision = replace(
            representative_relationship,
            relationship_type=RelationshipType.AMBIGUOUS_OVERLAP,
            status="ambiguous",
            reasons=incompatibilities,
        )
        representative: AlignmentHit | None = None
    else:
        decision = representative_relationship
        representative = representative_hit
    return PairRelationship(
        decision,
        representative,
        tuple(hit for hit, _ in accepted),
        tuple(rejected),
        incompatibilities,
    )


def _incompatibilities(
    accepted: list[tuple[AlignmentHit, Relationship]], tolerance: int
) -> tuple[str, ...]:
    first_hit, first_relationship = accepted[0]
    reasons: set[str] = set()
    for hit, relationship in accepted[1:]:
        if relationship.relationship_type != first_relationship.relationship_type:
            reasons.add(
                "accepted hits imply incompatible relationship types or terminal topologies"
            )
        if hit.orientation != first_hit.orientation:
            reasons.add("accepted hits imply incompatible orientations")
        placements = (
            abs(hit.query_start - first_hit.query_start),
            abs(hit.query_end - first_hit.query_end),
            abs(hit.target_start - first_hit.target_start),
            abs(hit.target_end - first_hit.target_end),
        )
        if any(distance > tolerance for distance in placements):
            reasons.add("accepted hits imply materially different coordinate placements")
    return tuple(sorted(reasons))


def _relationship(
    hit: AlignmentHit,
    relationship_type: RelationshipType,
    status: str,
    reasons: tuple[str, ...],
) -> Relationship:
    return Relationship(
        relationship_type=relationship_type,
        query_id=hit.query_id,
        target_id=hit.target_id,
        orientation=hit.orientation,
        identity=hit.identity,
        aligned_length=hit.alignment_block_length,
        query_coverage=hit.query_coverage,
        target_coverage=hit.target_coverage,
        status=status,
        reasons=reasons,
    )
