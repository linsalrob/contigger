"""Pseudomonas regression for metadata-only path planning."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from contigger.aligners.minimap2 import parse_paf
from contigger.config import build_run_config
from contigger.graph import build_relationship_graph
from contigger.models import (
    CatalogueMember,
    CatalogueSequence,
    Orientation,
    RelationshipGraph,
    RelationshipType,
    SequenceCatalogue,
)
from contigger.path_planning import plan_linear_paths
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.textio import open_text

BASELINE = Path(__file__).parents[1] / "benchmarks" / "pseudomonas_path_planning_baseline.json"
DATASET = Path(__file__).parents[1] / "test_data"


def graph_for_preset(preset: str) -> RelationshipGraph:
    """Build the source-ID diagnostic graph without cross-test imports."""
    with open_text(DATASET / "alignments" / f"all_vs_all.{preset}.paf.gz") as handle:
        decisions = tuple(
            classify_pair(group, build_run_config())
            for group in group_ordered_pairs(parse_paf(handle))
        )
    sequence_ids = tuple(
        sorted(
            {
                identifier
                for decision in decisions
                for identifier in (
                    decision.relationship.query_id,
                    decision.relationship.target_id,
                )
            }
        )
    )
    graph_decisions = tuple(
        decision
        for decision in decisions
        if decision.relationship.query_id != decision.relationship.target_id
        and decision.relationship.relationship_type is not RelationshipType.EXACT_MATCH
    )
    return build_relationship_graph(graph_decisions, sequence_ids)


@pytest.mark.parametrize("preset", ("asm5", "asm20"))
def test_pseudomonas_paths_remain_deferred_without_junction_support(preset: str) -> None:
    graph = graph_for_preset(preset)
    sequences = tuple(
        CatalogueSequence(node.sequence_id, "A", 1, node.sequence_id, node.sequence_id)
        for node in graph.nodes
    )
    members = tuple(
        CatalogueMember(
            node.sequence_id,
            node.sequence_id,
            "benchmark",
            node.sequence_id,
            Orientation.FORWARD,
            True,
        )
        for node in graph.nodes
    )
    result = plan_linear_paths(SequenceCatalogue(sequences, members), graph)
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))[preset]
    assert len(result.paths) == baseline["planned_paths"] == 0
    assert len(result.deferred_component_ids) == baseline["deferred_overlap_components"]
