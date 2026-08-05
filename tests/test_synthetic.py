"""Deterministic synthetic fixture and threshold boundary tests."""

from contigger.config import build_run_config
from contigger.models import AlignmentHit, Orientation, RelationshipType
from contigger.relationships import classify_alignment
from contigger.synthetic import mutate_substitutions, synthetic_cases


def test_synthetic_fixture_catalog_is_deterministic_and_complete() -> None:
    first = synthetic_cases()
    assert first == synthetic_cases()
    assert {case.name for case in first} == {
        "exact_duplicate",
        "reverse_complement_duplicate",
        "query_contained",
        "target_contained",
        "forward_suffix_prefix",
        "reverse_suffix_prefix",
        "identity_98_percent",
        "small_indel",
        "internal_conserved",
        "repeated_hits",
        "incompatible_terminal_placements",
        "low_complexity",
        "unrelated",
    }
    assert mutate_substitutions("AAAA", 2) == "CCAA"


def test_reverse_overlap_places_forward_target_match_at_suffix() -> None:
    case = next(case for case in synthetic_cases() if case.name == "reverse_suffix_prefix")
    reverse_complement_query_suffix = case.query[-300:]
    expected_forward_overlap = reverse_complement_query_suffix.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]
    assert case.target.endswith(expected_forward_overlap)


def _boundary_hit(*, end_distance: int = 10, matches: int = 98, block: int = 100) -> AlignmentHit:
    return AlignmentHit(
        query_id="q",
        target_id="t",
        query_length=210,
        target_length=210,
        query_start=210 - end_distance - block,
        query_end=210 - end_distance,
        target_start=end_distance,
        target_end=end_distance + block,
        orientation=Orientation.FORWARD,
        matching_bases=matches,
        alignment_block_length=block,
    )


def test_end_tolerance_boundary_is_inclusive() -> None:
    config = build_run_config(min_overlap=100, min_containment=100, end_tolerance=10)
    assert (
        classify_alignment(_boundary_hit(), config).relationship_type
        is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
    )
    assert (
        classify_alignment(_boundary_hit(end_distance=11), config).relationship_type
        is RelationshipType.NO_RELATIONSHIP
    )


def test_identity_and_length_threshold_boundaries() -> None:
    config = build_run_config(identity=98, min_overlap=100, min_containment=100, end_tolerance=10)
    assert (
        classify_alignment(_boundary_hit(), config).relationship_type
        is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
    )
    assert (
        classify_alignment(_boundary_hit(matches=97), config).relationship_type
        is RelationshipType.NO_RELATIONSHIP
    )
    short = _boundary_hit(matches=98, block=99)
    assert classify_alignment(short, config).relationship_type is RelationshipType.NO_RELATIONSHIP
