"""Checked-in Pseudomonas regression coverage for graph construction."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from contigger.aligners.minimap2 import parse_paf
from contigger.config import build_run_config
from contigger.graph import build_relationship_graph
from contigger.models import RelationshipGraph, RelationshipType
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.textio import open_text

DATASET = Path(__file__).parents[1] / "test_data"
BASELINE = DATASET.parent / "benchmarks" / "pseudomonas_graph_baseline.json"
AMBIGUOUS_PREFIXES = (
    "incompatible_placements_",
    "opposite_orientations_",
    "repeat_ambiguity_",
    "terminal_repeat_",
)


def graph_for_preset(preset: str) -> RelationshipGraph:
    """Build the unsimplified source-ID diagnostic graph for a checked-in PAF."""
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


@pytest.mark.parametrize(
    ("preset", "overlap_edges", "components"),
    (("asm5", 50, 58), ("asm20", 51, 57)),
)
def test_pseudomonas_graph_baseline(preset: str, overlap_edges: int, components: int) -> None:
    graph = graph_for_preset(preset)
    assert len(graph.nodes) == 90
    assert len(graph.containment_edges) == 3
    assert len(graph.overlap_edges) == overlap_edges
    assert graph.ambiguous_edges == ()
    assert len(graph.components) == components
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))[preset]
    assert baseline["nodes"] == len(graph.nodes)
    assert baseline["containment_edges"] == len(graph.containment_edges)
    assert baseline["overlap_edges"] == len(graph.overlap_edges)
    assert baseline["components"] == len(graph.components)
    ambiguous = [component for component in graph.components if component.ambiguous]
    assert len(ambiguous) == 1
    assert all(
        any(identifier.startswith(prefix) for identifier in ambiguous[0].sequence_ids)
        for prefix in AMBIGUOUS_PREFIXES
    )


def test_known_forbidden_boundary_remains_visible_but_is_not_a_merge() -> None:
    graph = graph_for_preset("asm20")
    pairs = {frozenset((edge.query_id, edge.target_id)) for edge in graph.overlap_edges}
    assert frozenset(("end_tolerance_51_left", "end_tolerance_51_right")) in pairs
