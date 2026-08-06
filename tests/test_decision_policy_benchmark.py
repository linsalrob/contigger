"""Checked-in Pseudomonas regression for conservative graph decisions."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from contigger.decision_policy import evaluate_graph_decisions
from contigger.models import GraphDecisionStatus
from tests.test_graph_benchmark import graph_for_preset

BASELINE = Path(__file__).parents[1] / "benchmarks" / "pseudomonas_decision_policy_baseline.json"


@pytest.mark.parametrize("preset", ("asm5", "asm20"))
def test_pseudomonas_decision_policy_baseline(preset: str) -> None:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))[preset]
    plan = evaluate_graph_decisions(graph_for_preset(preset))
    eligible_containments = tuple(
        item for item in plan.containment_decisions if item.status is GraphDecisionStatus.ELIGIBLE
    )
    deferred_containments = tuple(
        item for item in plan.containment_decisions if item.status is GraphDecisionStatus.DEFERRED
    )
    eligible_overlaps = tuple(
        item for item in plan.overlap_decisions if item.status is GraphDecisionStatus.ELIGIBLE
    )
    deferred_overlaps = tuple(
        item for item in plan.overlap_decisions if item.status is GraphDecisionStatus.DEFERRED
    )
    assert len(eligible_containments) == baseline["containment_eligible"]
    assert len(deferred_containments) == baseline["containment_deferred"]
    assert len(eligible_overlaps) == baseline["overlap_components_eligible"] == 0
    assert len(deferred_overlaps) == baseline["overlap_components_deferred"]
    assert (
        sum(len(item.edge_ids) for item in deferred_overlaps) == baseline["overlap_edges_deferred"]
    )


def test_known_forbidden_boundary_is_deferred_without_junction_support() -> None:
    graph = graph_for_preset("asm20")
    forbidden_edge = next(
        edge
        for edge in graph.overlap_edges
        if {edge.query_id, edge.target_id} == {"end_tolerance_51_left", "end_tolerance_51_right"}
    )
    plan = evaluate_graph_decisions(graph)
    decision = next(
        item for item in plan.overlap_decisions if forbidden_edge.edge_id in item.edge_ids
    )
    assert decision.status is GraphDecisionStatus.DEFERRED
    assert "lack explicit support" in " ".join(decision.reasons)
