"""Provenance-complete planning for eligible unambiguous linear paths."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from collections.abc import Iterable

from contigger.decision_policy import evaluate_graph_decisions
from contigger.exceptions import InputValidationError
from contigger.graph import validate_relationship_graph
from contigger.models import (
    CatalogueMember,
    GraphDecisionStatus,
    GraphEdge,
    LinearPathPlan,
    Orientation,
    PathPlanningResult,
    PlannedPathNode,
    PlannedSourceMember,
    RelationshipGraph,
    RelationshipType,
    SequenceCatalogue,
)


def plan_linear_paths(
    catalogue: SequenceCatalogue,
    graph: RelationshipGraph,
    *,
    junction_supported_edge_ids: Iterable[str] = (),
    intrinsically_safe_edge_ids: Iterable[str] = (),
) -> PathPlanningResult:
    """Plan canonical oriented paths without constructing or merging sequence."""
    members = _validate_catalogue(catalogue)
    validate_relationship_graph(graph)
    catalogue_ids = tuple(item.identifier for item in catalogue.sequences)
    graph_ids = tuple(item.sequence_id for item in graph.nodes)
    if graph_ids != catalogue_ids:
        missing = sorted(set(catalogue_ids) - set(graph_ids))
        unknown = sorted(set(graph_ids) - set(catalogue_ids))
        detail = f"missing {missing[0]}" if missing else f"unknown {unknown[0]}"
        raise InputValidationError(f"path-planning graph does not match catalogue: {detail}")

    decisions = evaluate_graph_decisions(
        graph,
        junction_supported_edge_ids=junction_supported_edge_ids,
        intrinsically_safe_edge_ids=intrinsically_safe_edge_ids,
    )
    edges = {edge.edge_id: edge for edge in graph.overlap_edges}
    paths = tuple(
        sorted(
            (
                _plan_component(decision.component_id, decision.edge_ids, edges, members)
                for decision in decisions.overlap_decisions
                if decision.status is GraphDecisionStatus.ELIGIBLE
            ),
            key=lambda item: item.path_id,
        )
    )
    deferred = tuple(
        decision.component_id
        for decision in decisions.overlap_decisions
        if decision.status is GraphDecisionStatus.DEFERRED
    )
    return PathPlanningResult(paths, tuple(sorted(deferred)))


def _validate_catalogue(
    catalogue: SequenceCatalogue,
) -> dict[str, tuple[CatalogueMember, ...]]:
    sequence_ids = tuple(item.identifier for item in catalogue.sequences)
    if sequence_ids != tuple(sorted(set(sequence_ids))):
        raise InputValidationError(
            "path planning requires uniquely and deterministically ordered catalogue sequences"
        )
    known = set(sequence_ids)
    grouped: dict[str, list[CatalogueMember]] = defaultdict(list)
    source_ids: set[str] = set()
    for member in catalogue.members:
        if member.catalogue_id not in known:
            raise InputValidationError(
                f"catalogue member references unknown sequence: {member.catalogue_id}"
            )
        if member.source_id in source_ids:
            raise InputValidationError(f"duplicate catalogue source member: {member.source_id}")
        source_ids.add(member.source_id)
        grouped[member.catalogue_id].append(member)
    missing = sorted(known - set(grouped))
    if missing:
        raise InputValidationError(f"catalogue sequence has no source provenance: {missing[0]}")
    return {
        sequence_id: tuple(sorted(items, key=lambda item: item.source_id))
        for sequence_id, items in grouped.items()
    }


def _plan_component(
    component_id: str,
    edge_ids: tuple[str, ...],
    edges: dict[str, GraphEdge],
    members: dict[str, tuple[CatalogueMember, ...]],
) -> LinearPathPlan:
    component_edges = tuple(edges[edge_id] for edge_id in edge_ids)
    adjacency: dict[str, list[GraphEdge]] = defaultdict(list)
    for edge in component_edges:
        adjacency[edge.query_id].append(edge)
        adjacency[edge.target_id].append(edge)
    endpoints = tuple(sorted(node for node, incident in adjacency.items() if len(incident) == 1))
    if len(endpoints) != 2 or any(len(incident) > 2 for incident in adjacency.values()):
        raise InputValidationError(f"eligible component is not a linear path: {component_id}")
    alternatives = tuple(_traverse(endpoint, adjacency, members) for endpoint in endpoints)
    nodes, ordered_edges = min(
        alternatives,
        key=lambda item: (
            tuple((node.sequence_id, node.orientation.value) for node in item[0]),
            item[1],
        ),
    )
    if len(ordered_edges) != len(component_edges):
        raise InputValidationError(f"eligible component path is disconnected: {component_id}")
    digest_input = "\0".join(
        [component_id]
        + [f"{node.sequence_id}:{node.orientation.value}" for node in nodes]
        + list(ordered_edges)
    )
    path_id = f"path_{hashlib.sha256(digest_input.encode()).hexdigest()[:20]}"
    return LinearPathPlan(path_id, component_id, nodes, ordered_edges)


def _traverse(
    start: str,
    adjacency: dict[str, list[GraphEdge]],
    members: dict[str, tuple[CatalogueMember, ...]],
) -> tuple[tuple[PlannedPathNode, ...], tuple[str, ...]]:
    nodes: list[PlannedPathNode] = []
    edge_ids: list[str] = []
    previous_edge: str | None = None
    current = start
    orientation: Orientation | None = None
    while True:
        remaining = [edge for edge in adjacency[current] if edge.edge_id != previous_edge]
        if not remaining:
            assert orientation is not None
            nodes.append(_planned_node(current, orientation, members[current]))
            break
        if len(remaining) != 1:
            raise InputValidationError(f"path traversal encountered a branch at {current}")
        edge = remaining[0]
        current_port, next_id, next_port = _edge_step(edge, current)
        if orientation is None:
            orientation = Orientation.FORWARD if current_port == "suffix" else Orientation.REVERSE
        elif _global_port(current_port, orientation) != "suffix":
            raise InputValidationError(f"path uses incompatible terminal at {current}")
        nodes.append(_planned_node(current, orientation, members[current]))
        next_orientation = Orientation.FORWARD if next_port == "prefix" else Orientation.REVERSE
        edge_ids.append(edge.edge_id)
        previous_edge = edge.edge_id
        current = next_id
        orientation = next_orientation
    return tuple(nodes), tuple(edge_ids)


def _planned_node(
    sequence_id: str,
    orientation: Orientation,
    members: tuple[CatalogueMember, ...],
) -> PlannedPathNode:
    planned_members = tuple(
        PlannedSourceMember(
            member.source_id,
            member.source_sample,
            member.original_identifier,
            _combine_orientation(member.orientation, orientation),
        )
        for member in members
    )
    return PlannedPathNode(sequence_id, orientation, planned_members)


def _combine_orientation(first: Orientation, second: Orientation) -> Orientation:
    return Orientation.FORWARD if first is second else Orientation.REVERSE


def _global_port(port: str, orientation: Orientation) -> str:
    if orientation is Orientation.FORWARD:
        return port
    return "suffix" if port == "prefix" else "prefix"


def _edge_step(edge: GraphEdge, current: str) -> tuple[str, str, str]:
    query_port, target_port = _edge_ports(edge)
    if current == edge.query_id:
        return query_port, edge.target_id, target_port
    return target_port, edge.query_id, query_port


def _edge_ports(edge: GraphEdge) -> tuple[str, str]:
    target_prefix = "prefix" if edge.orientation is Orientation.FORWARD else "suffix"
    target_suffix = "suffix" if edge.orientation is Orientation.FORWARD else "prefix"
    if edge.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX:
        return "suffix", target_prefix
    if edge.relationship_type is RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX:
        return "prefix", target_suffix
    raise InputValidationError(f"path edge is not a terminal overlap: {edge.edge_id}")
