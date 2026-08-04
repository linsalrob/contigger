"""Conservative pure-function relationship classification."""

from contigger.models import AlignmentHit, Relationship, RelationshipType, RunConfig


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
