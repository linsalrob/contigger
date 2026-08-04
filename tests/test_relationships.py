"""Conservative classification tests built without an aligner."""

import pytest

from contigger.config import build_run_config
from contigger.models import AlignmentHit, AlignmentType, Orientation, RelationshipType
from contigger.relationships import classify_alignment, classify_pair, group_ordered_pairs

CONFIG = build_run_config(min_overlap=100, min_containment=50, end_tolerance=10)


def hit(
    *,
    qlen: int = 1000,
    tlen: int = 1000,
    qs: int = 0,
    qe: int = 1000,
    ts: int = 0,
    te: int = 1000,
    matches: int | None = None,
    block: int | None = None,
    orientation: Orientation = Orientation.FORWARD,
) -> AlignmentHit:
    aligned = block if block is not None else max(qe - qs, te - ts)
    return AlignmentHit(
        query_id="q",
        target_id="t",
        query_length=qlen,
        target_length=tlen,
        query_start=qs,
        query_end=qe,
        target_start=ts,
        target_end=te,
        orientation=orientation,
        matching_bases=aligned if matches is None else matches,
        alignment_block_length=aligned,
    )


@pytest.mark.parametrize("orientation", [Orientation.FORWARD, Orientation.REVERSE])
def test_exact_match_in_both_orientations(orientation: Orientation) -> None:
    result = classify_alignment(hit(orientation=orientation), CONFIG)
    assert result.relationship_type is RelationshipType.EXACT_MATCH


def test_query_contained_in_target() -> None:
    result = classify_alignment(hit(qlen=500, tlen=1000, qe=500, ts=200, te=700, block=500), CONFIG)
    assert result.relationship_type is RelationshipType.QUERY_CONTAINED_IN_TARGET


def test_target_contained_in_query() -> None:
    result = classify_alignment(hit(qlen=1000, tlen=500, qs=200, qe=700, te=500, block=500), CONFIG)
    assert result.relationship_type is RelationshipType.TARGET_CONTAINED_IN_QUERY


def test_query_suffix_to_target_prefix() -> None:
    result = classify_alignment(hit(qs=700, qe=1000, ts=0, te=300, block=300), CONFIG)
    assert result.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX


def test_target_suffix_to_query_prefix() -> None:
    result = classify_alignment(hit(qs=0, qe=300, ts=700, te=1000, block=300), CONFIG)
    assert result.relationship_type is RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX


def test_high_identity_internal_alignment_is_not_relationship() -> None:
    result = classify_alignment(hit(qs=300, qe=600, ts=400, te=700, block=300), CONFIG)
    assert result.relationship_type is RelationshipType.NO_RELATIONSHIP
    assert "terminal geometry" in result.reasons[0]


def test_below_identity_threshold() -> None:
    result = classify_alignment(hit(qs=700, ts=0, te=300, block=300, matches=290), CONFIG)
    assert result.relationship_type is RelationshipType.NO_RELATIONSHIP
    assert "identity" in result.reasons[0]


def test_below_minimum_overlap() -> None:
    result = classify_alignment(hit(qs=950, qe=1000, ts=0, te=50, block=50), CONFIG)
    assert result.relationship_type is RelationshipType.NO_RELATIONSHIP
    assert "minimum overlap" in result.reasons[0]


def test_outside_end_tolerance() -> None:
    result = classify_alignment(hit(qs=700, qe=989, ts=11, te=300, block=289), CONFIG)
    assert result.relationship_type is RelationshipType.NO_RELATIONSHIP


def test_zero_length_alignment() -> None:
    result = classify_alignment(hit(qs=0, qe=0, ts=0, te=0, block=0, matches=0), CONFIG)
    assert result.relationship_type is RelationshipType.NO_RELATIONSHIP
    assert result.identity == 0.0


def test_ambiguous_competing_terminal_topology() -> None:
    broad_tolerance = build_run_config(min_overlap=50, min_containment=50, end_tolerance=60)
    result = classify_alignment(
        hit(qlen=200, tlen=200, qs=50, qe=150, ts=50, te=150, block=100), broad_tolerance
    )
    assert result.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP


def test_reverse_orientation_swaps_target_terminal_meaning() -> None:
    result = classify_alignment(
        hit(qs=700, qe=1000, ts=700, te=1000, block=300, orientation=Orientation.REVERSE),
        CONFIG,
    )
    assert result.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX


def test_duplicate_hits_are_ignored_deterministically() -> None:
    alignment = hit(qs=700, ts=0, te=300, block=300)
    result = classify_pair([alignment, alignment], CONFIG)
    assert result.relationship.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
    assert result.accepted_hits == (alignment,)


def test_equivalent_primary_and_secondary_hits_collapse() -> None:
    primary = hit(qs=700, ts=0, te=300, block=300)
    primary = AlignmentHit(
        **{
            field: getattr(primary, field)
            for field in primary.__dataclass_fields__
            if field != "alignment_type"
        },
        alignment_type=AlignmentType.PRIMARY,
    )
    secondary = AlignmentHit(
        **{
            field: getattr(primary, field)
            for field in primary.__dataclass_fields__
            if field
            not in {"query_start", "query_end", "target_start", "target_end", "alignment_type"}
        },
        query_start=705,
        query_end=995,
        target_start=5,
        target_end=295,
        alignment_type=AlignmentType.SECONDARY,
    )
    result = classify_pair([secondary, primary], CONFIG)
    assert result.relationship.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
    assert result.representative_hit == primary
    assert len(result.accepted_hits) == 2


def test_incompatible_secondary_hit_is_ambiguous() -> None:
    first = hit(qs=700, ts=0, te=300, block=300)
    second = hit(qs=0, qe=300, ts=700, te=1000, block=300)
    result = classify_pair([first, second], CONFIG)
    assert result.relationship.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP
    assert result.representative_hit is None
    assert any("topologies" in reason for reason in result.ambiguity_reasons)


def test_repeat_placements_are_ambiguous_even_with_same_topology() -> None:
    config = build_run_config(min_overlap=100, min_containment=50, end_tolerance=100)
    first = hit(qlen=1200, tlen=1200, qs=900, qe=1200, ts=0, te=300, block=300)
    second = hit(qlen=1200, tlen=1200, qs=750, qe=1100, ts=50, te=400, block=350)
    result = classify_pair([first, second], config)
    assert result.relationship.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP
    assert any("placements" in reason for reason in result.ambiguity_reasons)


def test_reverse_complement_orientation_conflict_is_ambiguous() -> None:
    forward = hit(qs=700, ts=0, te=300, block=300)
    reverse = hit(qs=700, ts=700, te=1000, block=300, orientation=Orientation.REVERSE)
    result = classify_pair([forward, reverse], CONFIG)
    assert result.relationship.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP
    assert any("orientations" in reason for reason in result.ambiguity_reasons)


def test_rejected_hits_and_diagnostics_are_retained() -> None:
    accepted = hit(qs=700, ts=0, te=300, block=300)
    internal = hit(qs=300, qe=600, ts=400, te=700, block=300)
    result = classify_pair([accepted, internal], CONFIG)
    assert result.relationship.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
    assert len(result.rejected_alignments) == 1
    assert "terminal geometry" in result.rejected_alignments[0].relationship.reasons[0]


def test_ordered_pair_grouping_is_deterministic() -> None:
    z_pair = hit()
    z_pair = AlignmentHit(
        **{
            field: getattr(z_pair, field)
            for field in z_pair.__dataclass_fields__
            if field not in {"query_id", "target_id"}
        },
        query_id="z",
        target_id="a",
    )
    a_pair = hit()
    groups = list(group_ordered_pairs([z_pair, a_pair]))
    assert [(group[0].query_id, group[0].target_id) for group in groups] == [("q", "t"), ("z", "a")]
