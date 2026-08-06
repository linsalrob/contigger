"""Tests for provenance-complete metadata-only linear path planning."""

from __future__ import annotations

from dataclasses import replace

import pytest

from contigger.exceptions import InputValidationError
from contigger.graph import build_relationship_graph
from contigger.models import (
    AlignmentHit,
    CatalogueMember,
    CatalogueSequence,
    Orientation,
    PairRelationship,
    Relationship,
    RelationshipType,
    SequenceCatalogue,
)
from contigger.path_planning import plan_linear_paths


def catalogue(*identifiers: str) -> SequenceCatalogue:
    """Build a small catalogue with two source members for the first sequence."""
    sequences = tuple(
        CatalogueSequence(identifier, "A" * 1000, 1000, identifier, f"S:{identifier}")
        for identifier in sorted(identifiers)
    )
    members = [
        CatalogueMember(identifier, f"S:{identifier}", "S", identifier, Orientation.FORWARD, True)
        for identifier in sorted(identifiers)
    ]
    members.append(
        CatalogueMember(identifiers[0], "T:duplicate", "T", "duplicate", Orientation.REVERSE)
    )
    return SequenceCatalogue(
        sequences, tuple(sorted(members, key=lambda item: (item.catalogue_id, item.source_id)))
    )


def overlap(
    query: str,
    target: str,
    *,
    orientation: Orientation = Orientation.FORWARD,
) -> PairRelationship:
    """Build one terminal-overlap decision."""
    target_start, target_end = (0, 200) if orientation is Orientation.FORWARD else (800, 1000)
    hit = AlignmentHit(
        query,
        target,
        1000,
        1000,
        800,
        1000,
        target_start,
        target_end,
        orientation,
        200,
        200,
    )
    relationship = Relationship(
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        query,
        target,
        orientation,
        1.0,
        200,
        0.2,
        0.2,
        "candidate",
    )
    return PairRelationship(relationship, hit, (hit,), ())


def test_supported_linear_path_has_canonical_order_and_complete_provenance() -> None:
    graph = build_relationship_graph((overlap("a", "b"), overlap("b", "c")), ("a", "b", "c"))
    supported = tuple(edge.edge_id for edge in graph.overlap_edges)
    plan = plan_linear_paths(
        catalogue("a", "b", "c"), graph, junction_supported_edge_ids=reversed(supported)
    )
    assert len(plan.paths) == 1
    path = plan.paths[0]
    assert tuple(node.sequence_id for node in path.nodes) == ("a", "b", "c")
    assert tuple(node.orientation for node in path.nodes) == (Orientation.FORWARD,) * 3
    assert path.edge_ids == supported
    assert {item.source_id for item in path.nodes[0].source_members} == {"S:a", "T:duplicate"}
    assert {item.orientation for item in path.nodes[0].source_members} == {
        Orientation.FORWARD,
        Orientation.REVERSE,
    }


def test_reverse_overlap_assigns_explicit_path_relative_orientation() -> None:
    graph = build_relationship_graph(
        (overlap("a", "b", orientation=Orientation.REVERSE),), ("a", "b")
    )
    edge_id = graph.overlap_edges[0].edge_id
    plan = plan_linear_paths(catalogue("a", "b"), graph, junction_supported_edge_ids=(edge_id,))
    orientations = {node.sequence_id: node.orientation for node in plan.paths[0].nodes}
    assert orientations["a"] is not orientations["b"]


def test_unsupported_and_ambiguous_components_remain_deferred() -> None:
    graph = build_relationship_graph((overlap("a", "b"), overlap("a", "c")), ("a", "b", "c"))
    plan = plan_linear_paths(catalogue("a", "b", "c"), graph)
    assert plan.paths == ()
    assert plan.deferred_component_ids == (graph.components[0].component_id,)


def test_planning_is_deterministic_across_relationship_order() -> None:
    decisions = (overlap("a", "b"), overlap("b", "c"))
    first_graph = build_relationship_graph(decisions, ("a", "b", "c"))
    second_graph = build_relationship_graph(reversed(decisions), ("c", "b", "a"))
    first = plan_linear_paths(
        catalogue("a", "b", "c"),
        first_graph,
        junction_supported_edge_ids=(edge.edge_id for edge in first_graph.overlap_edges),
    )
    second = plan_linear_paths(
        catalogue("a", "b", "c"),
        second_graph,
        junction_supported_edge_ids=(edge.edge_id for edge in second_graph.overlap_edges),
    )
    assert first == second


def test_catalogue_graph_mismatch_and_missing_provenance_are_rejected() -> None:
    graph = build_relationship_graph((overlap("a", "b"),), ("a", "b"))
    with pytest.raises(InputValidationError, match="does not match catalogue: missing c"):
        plan_linear_paths(catalogue("a", "b", "c"), graph)
    incomplete = replace(
        catalogue("a", "b"),
        members=tuple(item for item in catalogue("a", "b").members if item.catalogue_id != "b"),
    )
    with pytest.raises(InputValidationError, match="no source provenance: b"):
        plan_linear_paths(incomplete, graph)
