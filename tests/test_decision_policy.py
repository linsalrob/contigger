"""Tests for conservative graph decision policy."""

from __future__ import annotations

from dataclasses import replace

import pytest

from contigger.decision_policy import evaluate_graph_decisions
from contigger.exceptions import InputValidationError
from contigger.graph import build_relationship_graph
from contigger.models import GraphDecisionStatus, RelationshipType
from tests.test_graph import decision


def test_unique_containment_is_eligible_without_removing_a_node() -> None:
    graph = build_relationship_graph(
        (
            decision(
                "short",
                "long",
                RelationshipType.QUERY_CONTAINED_IN_TARGET,
                query_start=0,
                query_end=1000,
                target_start=0,
                target_end=1000,
            ),
        )
    )
    plan = evaluate_graph_decisions(graph)
    assert plan.containment_decisions[0].status is GraphDecisionStatus.ELIGIBLE
    assert plan.containment_decisions[0].contained_sequence_id == "short"
    assert [node.sequence_id for node in graph.nodes] == ["long", "short"]


def test_multiple_containers_are_deferred() -> None:
    graph = build_relationship_graph(
        (
            decision("short", "long_a", RelationshipType.QUERY_CONTAINED_IN_TARGET),
            decision("short", "long_b", RelationshipType.QUERY_CONTAINED_IN_TARGET),
        )
    )
    plan = evaluate_graph_decisions(graph)
    assert {item.status for item in plan.containment_decisions} == {GraphDecisionStatus.DEFERRED}
    assert all(
        "multiple possible containers" in " ".join(item.reasons)
        for item in plan.containment_decisions
    )


def test_overlap_without_junction_support_is_deferred() -> None:
    graph = build_relationship_graph(
        (decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),)
    )
    plan = evaluate_graph_decisions(graph)
    assert plan.overlap_decisions[0].status is GraphDecisionStatus.DEFERRED
    assert (
        "one or more proposed junctions lack explicit support" in plan.overlap_decisions[0].reasons
    )


def test_supported_unambiguous_overlap_component_is_eligible_for_later_planning() -> None:
    graph = build_relationship_graph(
        (
            decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("b", "c", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
        )
    )
    edge_ids = tuple(edge.edge_id for edge in graph.overlap_edges)
    plan = evaluate_graph_decisions(graph, junction_supported_edge_ids=reversed(edge_ids))
    assert plan.overlap_decisions[0].status is GraphDecisionStatus.ELIGIBLE
    assert plan.overlap_decisions[0].edge_ids == tuple(sorted(edge_ids))


def test_supported_branch_remains_deferred() -> None:
    graph = build_relationship_graph(
        (
            decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
            decision("a", "c", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),
        )
    )
    plan = evaluate_graph_decisions(
        graph,
        junction_supported_edge_ids=(edge.edge_id for edge in graph.overlap_edges),
    )
    assert plan.overlap_decisions[0].status is GraphDecisionStatus.DEFERRED
    assert "multiple overlaps compete" in " ".join(plan.overlap_decisions[0].reasons)


def test_unknown_or_duplicate_supported_edges_are_rejected() -> None:
    graph = build_relationship_graph(
        (decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),)
    )
    edge_id = graph.overlap_edges[0].edge_id
    with pytest.raises(InputValidationError, match="must be unique"):
        evaluate_graph_decisions(graph, junction_supported_edge_ids=(edge_id, edge_id))
    with pytest.raises(InputValidationError, match="unknown overlap edge"):
        evaluate_graph_decisions(graph, junction_supported_edge_ids=("missing",))


def test_inconsistent_graph_component_membership_is_rejected() -> None:
    graph = build_relationship_graph(
        (decision("a", "b", RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX),)
    )
    broken_component = replace(graph.components[0], sequence_ids=("a",))
    broken_graph = replace(graph, components=(broken_component,))
    with pytest.raises(InputValidationError, match="absent from components: b"):
        evaluate_graph_decisions(broken_graph)
