"""Conservative sequence-construction regression tests."""

import pytest

from contigger.exceptions import InputValidationError
from contigger.merge import (
    _construct_path,
    _edge_has_exact_reconcilable_overlap,
    _relationships_for_graph,
)
from contigger.models import (
    AlignmentHit,
    CatalogueSequence,
    GraphEdge,
    GraphEdgeKind,
    LinearPathPlan,
    Orientation,
    PairRelationship,
    PlannedPathNode,
    Relationship,
    RelationshipType,
    SequenceCatalogue,
)


def _path(edge: GraphEdge) -> LinearPathPlan:
    """Build a two-node forward path for a synthetic exact overlap."""
    return LinearPathPlan(
        "path_test",
        "component_test",
        (
            PlannedPathNode("a", Orientation.FORWARD, ()),
            PlannedPathNode("b", Orientation.FORWARD, ()),
        ),
        (edge.edge_id,),
    )


def test_constructs_perfect_forward_terminal_overlap() -> None:
    edge = GraphEdge(
        "edge_ab",
        GraphEdgeKind.OVERLAP,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        "a",
        "b",
        Orientation.FORWARD,
        4,
        8,
        0,
        4,
        1.0,
        4,
        1,
        0,
    )
    sequence, intervals = _construct_path(
        _path(edge),
        {edge.edge_id: edge},
        {
            "a": CatalogueSequence("a", "AAAACCCC", 8, "a", "S:a"),
            "b": CatalogueSequence("b", "CCCCGGGG", 8, "b", "S:b"),
        },
    )
    assert sequence == "AAAACCCCGGGG"
    assert intervals == ((0, 8), (4, 12))


def test_constructs_perfect_reverse_complement_terminal_overlap() -> None:
    edge = GraphEdge(
        "edge_ab_rc",
        GraphEdgeKind.OVERLAP,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        "a",
        "b",
        Orientation.REVERSE,
        4,
        8,
        4,
        8,
        1.0,
        4,
        1,
        0,
    )
    path = LinearPathPlan(
        "path_rc",
        "component_rc",
        (
            PlannedPathNode("a", Orientation.FORWARD, ()),
            PlannedPathNode("b", Orientation.REVERSE, ()),
        ),
        (edge.edge_id,),
    )
    sequence, _ = _construct_path(
        path,
        {edge.edge_id: edge},
        {
            "a": CatalogueSequence("a", "AAAAGGGG", 8, "a", "S:a"),
            "b": CatalogueSequence("b", "AAAACCCC", 8, "b", "S:b"),
        },
    )
    assert sequence == "AAAAGGGGTTTT"


def test_imperfect_overlap_is_rejected_without_consensus() -> None:
    edge = GraphEdge(
        "edge_ab",
        GraphEdgeKind.OVERLAP,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        "a",
        "b",
        Orientation.FORWARD,
        4,
        8,
        0,
        4,
        0.75,
        4,
        1,
        0,
    )
    with pytest.raises(InputValidationError, match="imperfect"):
        _construct_path(
            _path(edge),
            {edge.edge_id: edge},
            {
                "a": CatalogueSequence("a", "AAAACCCC", 8, "a", "S:a"),
                "b": CatalogueSequence("b", "CCCTGGGG", 8, "b", "S:b"),
            },
        )


def test_terminal_overhang_is_deferred_instead_of_overwriting_the_junction() -> None:
    edge = GraphEdge(
        "edge_overhang",
        GraphEdgeKind.OVERLAP,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        "a",
        "b",
        Orientation.FORWARD,
        4,
        8,
        0,
        4,
        1.0,
        4,
        1,
        0,
    )
    with pytest.raises(InputValidationError, match="terminal"):
        _construct_path(
            _path(edge),
            {edge.edge_id: edge},
            {
                "a": CatalogueSequence("a", "AAAACCCCX", 9, "a", "S:a"),
                "b": CatalogueSequence("b", "CCCCGGGG", 8, "b", "S:b"),
            },
        )


def test_endpoint_beyond_tolerance_cannot_become_a_known_false_join() -> None:
    edge = GraphEdge(
        "forbidden",
        GraphEdgeKind.OVERLAP,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        "a",
        "b",
        Orientation.FORWARD,
        0,
        7,
        0,
        4,
        1.0,
        4,
        1,
        0,
    )
    catalogue = SequenceCatalogue(
        (
            CatalogueSequence("a", "AAAACCCC", 8, "a", "S:a"),
            CatalogueSequence("b", "CCCCGGGG", 8, "b", "S:b"),
        ),
        (),
    )
    assert not _edge_has_exact_reconcilable_overlap(edge, catalogue, 0)


def test_exact_paf_relationship_is_retained_only_outside_the_graph() -> None:
    """Exact PAF classifications cannot bypass catalogue-level deduplication."""
    hit = AlignmentHit("a", "b", 8, 8, 0, 8, 0, 8, Orientation.FORWARD, 8, 8)
    exact = PairRelationship(
        Relationship(
            RelationshipType.EXACT_MATCH,
            "a",
            "b",
            Orientation.FORWARD,
            1.0,
            8,
            1.0,
            1.0,
            "confident",
            ("full exact span",),
        ),
        hit,
        (hit,),
        (),
    )
    assert _relationships_for_graph((exact,)) == ()
