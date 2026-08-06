"""Ambiguity-preserving relationship graph tests."""

from __future__ import annotations

import pytest

from contigger.exceptions import InputValidationError
from contigger.graph import build_components, build_relationship_graph, validate_relationship_graph
from contigger.models import (
    AlignmentHit,
    GraphEdgeKind,
    Orientation,
    PairRelationship,
    RejectedAlignment,
    Relationship,
    RelationshipType,
)


def decision(
    query: str,
    target: str,
    relationship_type: RelationshipType,
    *,
    orientation: Orientation = Orientation.FORWARD,
    query_start: int = 800,
    query_end: int = 1000,
    target_start: int = 0,
    target_end: int = 200,
) -> PairRelationship:
    """Build a complete synthetic pair decision with representative coordinates."""
    ambiguous = relationship_type is RelationshipType.AMBIGUOUS_OVERLAP
    rejected = relationship_type is RelationshipType.NO_RELATIONSHIP
    query_span = query_end - query_start
    target_span = target_end - target_start
    hit = AlignmentHit(
        query,
        target,
        1000,
        1000,
        query_start,
        query_end,
        target_start,
        target_end,
        orientation,
        min(query_span, target_span),
        max(query_span, target_span),
    )
    relationship = Relationship(
        relationship_type,
        query,
        target,
        orientation,
        1.0,
        200,
        0.2,
        0.2,
        "ambiguous" if ambiguous else "rejected" if rejected else "candidate",
        ("synthetic decision",),
    )
    return PairRelationship(
        relationship,
        None if ambiguous else hit,
        (hit,) if not rejected else (),
        (),
        ("competing placements",) if ambiguous else (),
    )


def reciprocal(item: PairRelationship) -> PairRelationship:
    """Return the geometrically equivalent inverse ordered decision."""
    inverse = {
        RelationshipType.QUERY_CONTAINED_IN_TARGET: RelationshipType.TARGET_CONTAINED_IN_QUERY,
        RelationshipType.TARGET_CONTAINED_IN_QUERY: RelationshipType.QUERY_CONTAINED_IN_TARGET,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX: (
            RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX
        ),
        RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX: (
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
        ),
    }
    hit = item.representative_hit
    assert hit is not None
    relationship_type = inverse[item.relationship.relationship_type]
    if hit.orientation is Orientation.REVERSE and item.relationship.relationship_type in {
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX,
    }:
        relationship_type = item.relationship.relationship_type
    return decision(
        hit.target_id,
        hit.query_id,
        relationship_type,
        orientation=hit.orientation,
        query_start=hit.target_start,
        query_end=hit.target_end,
        target_start=hit.query_start,
        target_end=hit.query_end,
    )


def test_reciprocal_overlap_collapses_to_one_deterministic_edge() -> None:
    forward = decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX)
    graph = build_relationship_graph((reciprocal(forward), forward), ("b", "a", "isolated"))
    assert [node.sequence_id for node in graph.nodes] == ["a", "b", "isolated"]
    assert len(graph.overlap_edges) == 1
    edge = graph.overlap_edges[0]
    assert edge.kind is GraphEdgeKind.OVERLAP
    assert (edge.query_id, edge.target_id) == ("a", "b")
    assert edge.query_start == 800
    assert edge.accepted_hit_count == 2
    assert edge.rejected_hit_count == 0
    assert graph.containment_edges == graph.ambiguous_edges == ()
    assert [component.sequence_ids for component in graph.components] == [
        ("a", "b"),
        ("isolated",),
    ]
    assert not graph.components[0].ambiguous
    validate_relationship_graph(graph)


def test_containment_is_structurally_separate_from_overlap() -> None:
    contained = decision(
        "short",
        "long",
        RelationshipType.QUERY_CONTAINED_IN_TARGET,
        query_start=0,
        query_end=1000,
        target_start=0,
        target_end=1000,
    )
    graph = build_relationship_graph((contained, reciprocal(contained)))
    assert len(graph.containment_edges) == 1
    assert graph.overlap_edges == ()
    assert not graph.components[0].ambiguous


def test_reverse_complement_reciprocals_preserve_query_relative_topology() -> None:
    forward = decision(
        "a",
        "b",
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        orientation=Orientation.REVERSE,
        target_start=800,
        target_end=1000,
    )
    graph = build_relationship_graph((forward, reciprocal(forward)))
    assert len(graph.overlap_edges) == 1
    assert graph.ambiguous_edges == ()
    assert graph.overlap_edges[0].orientation is Orientation.REVERSE


def test_no_relationship_is_not_an_edge_but_nodes_remain() -> None:
    rejected = decision("a", "b", RelationshipType.NO_RELATIONSHIP)
    graph = build_relationship_graph((rejected,), ("a", "b"))
    assert graph.containment_edges == graph.overlap_edges == graph.ambiguous_edges == ()
    assert [component.sequence_ids for component in graph.components] == [("a",), ("b",)]


def test_pair_ambiguity_is_retained_as_an_ambiguous_edge() -> None:
    ambiguous = decision("a", "b", RelationshipType.AMBIGUOUS_OVERLAP)
    graph = build_relationship_graph((ambiguous,))
    assert len(graph.ambiguous_edges) == 1
    assert graph.ambiguous_edges[0].query_start is None
    assert graph.components[0].ambiguous
    assert "ambiguous pairwise evidence" in graph.components[0].ambiguity_reasons[0]


def test_terminal_branch_marks_component_ambiguous_without_dropping_edges() -> None:
    graph = build_relationship_graph(
        (
            decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("a", "c", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
        )
    )
    assert len(graph.overlap_edges) == 2
    assert graph.components[0].ambiguous
    assert (
        "multiple overlaps compete for the same oriented terminal"
        in graph.components[0].ambiguity_reasons
    )


def test_degree_two_linear_path_is_not_mistaken_for_a_branch() -> None:
    graph = build_relationship_graph(
        (
            decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("b", "c", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
        )
    )
    assert not graph.components[0].ambiguous


def test_multiple_containers_are_preserved_and_marked_ambiguous() -> None:
    graph = build_relationship_graph(
        (
            decision("short", "long_a", RelationshipType.QUERY_CONTAINED_IN_TARGET),
            decision("short", "long_b", RelationshipType.QUERY_CONTAINED_IN_TARGET),
        )
    )
    assert len(graph.containment_edges) == 2
    assert "multiple possible containers" in graph.components[0].ambiguity_reasons[0]


def test_inconsistent_reciprocals_become_ambiguous_instead_of_electing_a_winner() -> None:
    first = decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX)
    conflicting = decision(
        "b",
        "a",
        RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX,
        query_start=0,
        query_end=200,
        target_start=799,
        target_end=999,
    )
    graph = build_relationship_graph((first, conflicting))
    assert graph.overlap_edges == ()
    assert len(graph.ambiguous_edges) == 1
    assert "inconsistent" in " ".join(graph.ambiguous_edges[0].reasons)


def test_relationship_mismatch_retains_rejected_reciprocal_diagnostics() -> None:
    overlap = decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX)
    rejected_base = decision("b", "a", RelationshipType.NO_RELATIONSHIP)
    assert rejected_base.representative_hit is not None
    rejected_relationship = Relationship(
        RelationshipType.NO_RELATIONSHIP,
        "b",
        "a",
        Orientation.FORWARD,
        1.0,
        200,
        0.2,
        0.2,
        "rejected",
        ("rejected reciprocal evidence",),
    )
    rejected = PairRelationship(
        rejected_relationship,
        rejected_base.representative_hit,
        (),
        (RejectedAlignment(rejected_base.representative_hit, rejected_relationship),),
    )
    graph = build_relationship_graph((overlap, rejected))
    edge = graph.ambiguous_edges[0]
    assert edge.accepted_hit_count == 1
    assert edge.rejected_hit_count == 1
    assert "rejected reciprocal evidence" in edge.reasons


def test_cycle_and_orientation_conflict_are_reported_without_simplification() -> None:
    graph = build_relationship_graph(
        (
            decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("b", "c", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision(
                "c",
                "a",
                RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
                orientation=Orientation.REVERSE,
            ),
        )
    )
    assert len(graph.overlap_edges) == 3
    assert "overlap subgraph contains a cycle" in graph.components[0].ambiguity_reasons
    assert "overlap orientations are mutually inconsistent" in graph.components[0].ambiguity_reasons


def test_overlap_cycle_detection_handles_disconnected_overlap_subgraphs() -> None:
    graph = build_relationship_graph(
        (
            decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("b", "c", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("c", "a", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("d", "e", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision(
                "c",
                "d",
                RelationshipType.QUERY_CONTAINED_IN_TARGET,
                query_start=0,
                query_end=1000,
                target_start=0,
                target_end=1000,
            ),
        )
    )
    assert "overlap subgraph contains a cycle" in graph.components[0].ambiguity_reasons


def test_invalid_graph_inputs_fail_clearly() -> None:
    exact = decision("a", "b", RelationshipType.EXACT_MATCH)
    with pytest.raises(InputValidationError, match="sequence catalogue"):
        build_relationship_graph((exact,))
    relationship = decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX)
    with pytest.raises(InputValidationError, match="duplicate ordered"):
        build_relationship_graph((relationship, relationship))
    with pytest.raises(InputValidationError, match="unknown sequence"):
        build_relationship_graph((relationship,), ("a",))
    mismatched_orientation = PairRelationship(
        Relationship(
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
            "a",
            "b",
            Orientation.REVERSE,
            1.0,
            200,
            0.2,
            0.2,
            "candidate",
            ("synthetic decision",),
        ),
        relationship.representative_hit,
        relationship.accepted_hits,
        relationship.rejected_alignments,
        (),
    )
    with pytest.raises(InputValidationError, match="orientation"):
        build_relationship_graph((mismatched_orientation,))


def test_component_wrapper_and_graph_order_are_deterministic() -> None:
    relationships = (
        decision("z", "y", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
        decision("b", "a", RelationshipType.QUERY_CONTAINED_IN_TARGET),
    )
    first = build_relationship_graph(relationships)
    second = build_relationship_graph(tuple(reversed(relationships)))
    assert first == second
    assert build_components(relationships) == first.components
