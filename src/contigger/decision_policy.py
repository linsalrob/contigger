"""Conservative policy decisions over an unsimplified relationship graph."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable

from contigger.exceptions import InputValidationError
from contigger.graph import validate_relationship_graph
from contigger.models import (
    ContainmentDecision,
    GraphDecisionPlan,
    GraphDecisionStatus,
    GraphEdge,
    GraphEdgeKind,
    MergeComponent,
    OverlapComponentDecision,
    RelationshipGraph,
    RelationshipType,
)


def evaluate_graph_decisions(
    graph: RelationshipGraph,
    *,
    junction_supported_edge_ids: Iterable[str] = (),
) -> GraphDecisionPlan:
    """Evaluate conservative containment and overlap eligibility.

    Eligibility is only a typed input to later provenance-complete planning. It
    never removes a node, selects sequence, or authorizes a biological merge.
    Every overlap edge requires explicit junction support; graph structure and
    alignment identity alone are insufficient.
    """
    validate_relationship_graph(graph)
    components, component_by_node = _component_lookups(graph)
    overlap_ids = {edge.edge_id for edge in graph.overlap_edges}
    supported = tuple(sorted(junction_supported_edge_ids))
    if len(set(supported)) != len(supported):
        raise InputValidationError("junction-supported graph edge identifiers must be unique")
    unknown = sorted(set(supported) - overlap_ids)
    if unknown:
        raise InputValidationError(
            f"junction support references unknown overlap edge identifier: {unknown[0]}"
        )

    containment_by_child: dict[str, list[GraphEdge]] = defaultdict(list)
    for edge in graph.containment_edges:
        contained, _ = _contained_and_container(edge)
        containment_by_child[contained].append(edge)

    containment_decisions = tuple(
        _containment_decision(
            edge,
            components[component_by_node[_contained_and_container(edge)[0]]],
            containment_by_child,
        )
        for edge in sorted(graph.containment_edges, key=lambda item: item.edge_id)
    )
    supported_set = set(supported)
    overlap_decisions = tuple(
        decision
        for component in sorted(graph.components, key=lambda item: item.sequence_ids)
        if (decision := _overlap_decision(component, graph, supported_set)) is not None
    )
    return GraphDecisionPlan(containment_decisions, overlap_decisions)


def _component_lookups(
    graph: RelationshipGraph,
) -> tuple[dict[str, MergeComponent], dict[str, str]]:
    """Build lookups after public structural graph validation."""
    components: dict[str, MergeComponent] = {}
    component_by_node: dict[str, str] = {}
    for component in graph.components:
        components[component.component_id] = component
        for sequence_id in component.sequence_ids:
            component_by_node[sequence_id] = component.component_id
    return components, component_by_node


def _containment_decision(
    edge: GraphEdge,
    component: MergeComponent,
    containment_by_child: dict[str, list[GraphEdge]],
) -> ContainmentDecision:
    contained, container = _contained_and_container(edge)
    reasons: set[str] = set()
    if component.ambiguous:
        reasons.add("component ambiguity prevents containment disposition")
        reasons.update(component.ambiguity_reasons)
    if len(containment_by_child[contained]) != 1:
        reasons.add("contained sequence has multiple possible containers")
    status = GraphDecisionStatus.DEFERRED if reasons else GraphDecisionStatus.ELIGIBLE
    if not reasons:
        reasons.add("unique containment creates no new sequence junction")
    return ContainmentDecision(
        edge.edge_id,
        contained,
        container,
        status,
        tuple(sorted(reasons)),
    )


def _overlap_decision(
    component: MergeComponent,
    graph: RelationshipGraph,
    supported: set[str],
) -> OverlapComponentDecision | None:
    component_edge_ids = set(component.relationship_ids)
    overlap_ids = tuple(
        sorted(edge.edge_id for edge in graph.overlap_edges if edge.edge_id in component_edge_ids)
    )
    if not overlap_ids:
        return None
    reasons: set[str] = set()
    if component.ambiguous:
        reasons.add("component ambiguity prevents overlap path eligibility")
        reasons.update(component.ambiguity_reasons)
    if any(edge.edge_id in component_edge_ids for edge in graph.containment_edges):
        reasons.add("containment and overlap evidence require provenance-complete joint planning")
    unsupported = tuple(edge_id for edge_id in overlap_ids if edge_id not in supported)
    if unsupported:
        reasons.add("one or more proposed junctions lack explicit support")
    status = GraphDecisionStatus.DEFERRED if reasons else GraphDecisionStatus.ELIGIBLE
    if not reasons:
        reasons.add("all junctions are supported in an unambiguous overlap-only component")
    return OverlapComponentDecision(
        component.component_id,
        component.sequence_ids,
        overlap_ids,
        status,
        tuple(sorted(reasons)),
    )


def _contained_and_container(edge: GraphEdge) -> tuple[str, str]:
    if edge.kind is not GraphEdgeKind.CONTAINMENT:
        raise InputValidationError(f"edge is not a containment relationship: {edge.edge_id}")
    if edge.relationship_type is RelationshipType.QUERY_CONTAINED_IN_TARGET:
        return edge.query_id, edge.target_id
    if edge.relationship_type is RelationshipType.TARGET_CONTAINED_IN_QUERY:
        return edge.target_id, edge.query_id
    raise InputValidationError(
        f"containment edge has incompatible relationship type: {edge.edge_id}"
    )
