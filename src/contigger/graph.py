"""Deterministic ambiguity-preserving relationship graph construction."""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict, deque
from collections.abc import Iterable

from contigger.exceptions import InputValidationError
from contigger.models import (
    GraphEdge,
    GraphEdgeKind,
    GraphNode,
    MergeComponent,
    Orientation,
    PairRelationship,
    RelationshipGraph,
    RelationshipType,
)

_CONTAINMENTS = frozenset(
    {
        RelationshipType.QUERY_CONTAINED_IN_TARGET,
        RelationshipType.TARGET_CONTAINED_IN_QUERY,
    }
)
_OVERLAPS = frozenset(
    {
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX,
    }
)


def build_relationship_graph(
    relationships: Iterable[PairRelationship],
    sequence_ids: Iterable[str] = (),
) -> RelationshipGraph:
    """Build an unsimplified graph from complete ordered-pair decisions.

    Reciprocal observations collapse only when their types, orientations, and
    representative coordinates agree exactly after swapping query and target.
    Conflicts become explicit ambiguous edges; no score elects a winner.
    """
    decisions = tuple(
        sorted(
            relationships,
            key=lambda item: (
                item.relationship.query_id,
                item.relationship.target_id,
            ),
        )
    )
    supplied_nodes = tuple(sorted(sequence_ids))
    if any(not identifier for identifier in supplied_nodes):
        raise InputValidationError("graph sequence identifiers cannot be empty")
    if len(set(supplied_nodes)) != len(supplied_nodes):
        raise InputValidationError("graph sequence identifiers must be unique")

    ordered_pairs: set[tuple[str, str]] = set()
    observed_nodes: set[str] = set()
    grouped: dict[tuple[str, str], list[PairRelationship]] = defaultdict(list)
    for decision in decisions:
        relationship = decision.relationship
        pair = (relationship.query_id, relationship.target_id)
        if not pair[0] or not pair[1]:
            raise InputValidationError("graph relationships require non-empty identifiers")
        if pair[0] == pair[1]:
            raise InputValidationError(f"graph relationship cannot be a self-pair: {pair[0]}")
        if pair in ordered_pairs:
            raise InputValidationError(f"duplicate ordered graph relationship: {pair}")
        ordered_pairs.add(pair)
        observed_nodes.update(pair)
        if relationship.relationship_type is RelationshipType.EXACT_MATCH:
            raise InputValidationError(
                "exact matches must be resolved by the sequence catalogue before graph construction"
            )
        _validate_representative(decision)
        first, second = sorted(pair)
        grouped[(first, second)].append(decision)

    node_ids = tuple(sorted(set(supplied_nodes) | observed_nodes))
    if supplied_nodes:
        unknown = sorted(observed_nodes - set(supplied_nodes))
        if unknown:
            raise InputValidationError(
                f"graph relationship references unknown sequence identifier: {unknown[0]}"
            )

    edges = tuple(
        edge
        for pair in sorted(grouped)
        if (edge := _canonical_edge(pair, grouped[pair])) is not None
    )
    containment = tuple(edge for edge in edges if edge.kind is GraphEdgeKind.CONTAINMENT)
    overlap = tuple(edge for edge in edges if edge.kind is GraphEdgeKind.OVERLAP)
    ambiguous = tuple(edge for edge in edges if edge.kind is GraphEdgeKind.AMBIGUOUS)
    components = _components(node_ids, edges)
    return RelationshipGraph(
        nodes=tuple(GraphNode(identifier) for identifier in node_ids),
        containment_edges=containment,
        overlap_edges=overlap,
        ambiguous_edges=ambiguous,
        components=components,
    )


def build_components(
    relationships: tuple[PairRelationship, ...],
) -> tuple[MergeComponent, ...]:
    """Return deterministic components without simplifying or merging them."""
    return build_relationship_graph(relationships).components


def _validate_representative(decision: PairRelationship) -> None:
    relationship = decision.relationship
    hit = decision.representative_hit
    merge_like = relationship.relationship_type in _CONTAINMENTS | _OVERLAPS
    if merge_like and hit is None:
        raise InputValidationError(
            "non-ambiguous graph relationship requires a representative alignment: "
            f"{relationship.query_id}, {relationship.target_id}"
        )
    if hit is not None and (hit.query_id, hit.target_id) != (
        relationship.query_id,
        relationship.target_id,
    ):
        raise InputValidationError(
            "graph representative alignment identifiers do not match its relationship"
        )
    if hit is not None and hit.orientation is not relationship.orientation:
        raise InputValidationError(
            "graph representative alignment orientation does not match its relationship"
        )


def _canonical_edge(pair: tuple[str, str], decisions: list[PairRelationship]) -> GraphEdge | None:
    active = [
        decision
        for decision in decisions
        if decision.relationship.relationship_type is not RelationshipType.NO_RELATIONSHIP
    ]
    if not active:
        return None
    reasons: set[str] = set()
    if len(active) != len(decisions):
        reasons.add("reciprocal decisions disagree on whether a relationship exists")
    if any(
        decision.relationship.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP
        for decision in active
    ):
        reasons.add("pair classifier retained ambiguous alignment evidence")
    if len(active) == 2 and not _reciprocals_agree(active[0], active[1]):
        reasons.add("reciprocal ordered-pair decisions are inconsistent")
    if len(active) > 2:
        raise InputValidationError(f"too many ordered graph relationships for pair: {pair}")

    selected = min(
        active,
        key=lambda item: (item.relationship.query_id, item.relationship.target_id),
    )
    if reasons:
        return _edge_from_decision(
            pair,
            selected,
            GraphEdgeKind.AMBIGUOUS,
            RelationshipType.AMBIGUOUS_OVERLAP,
            tuple(sorted(reasons | _decision_reasons(decisions))),
            retain_coordinates=False,
            accepted_hit_count=sum(len(item.accepted_hits) for item in decisions),
            rejected_hit_count=sum(len(item.rejected_alignments) for item in decisions),
        )

    relationship_type = selected.relationship.relationship_type
    if relationship_type in _CONTAINMENTS:
        kind = GraphEdgeKind.CONTAINMENT
    elif relationship_type in _OVERLAPS:
        kind = GraphEdgeKind.OVERLAP
    else:
        raise InputValidationError(f"unsupported graph relationship type: {relationship_type}")
    return _edge_from_decision(
        pair,
        selected,
        kind,
        relationship_type,
        tuple(sorted(_decision_reasons(active))),
        retain_coordinates=True,
        accepted_hit_count=sum(len(item.accepted_hits) for item in active),
        rejected_hit_count=sum(len(item.rejected_alignments) for item in active),
    )


def _reciprocals_agree(first: PairRelationship, second: PairRelationship) -> bool:
    left = first.relationship
    right = second.relationship
    if (left.query_id, left.target_id) != (right.target_id, right.query_id):
        return False
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
    reciprocal_type = inverse.get(left.relationship_type)
    if left.orientation is Orientation.REVERSE and left.relationship_type in _OVERLAPS:
        # Reversing a reverse-strand ordered pair preserves the physical terminal
        # on both sequences, so the named query-relative topology also stays the same.
        reciprocal_type = left.relationship_type
    if reciprocal_type is not right.relationship_type:
        return False
    if left.orientation is not right.orientation:
        return False
    first_hit = first.representative_hit
    second_hit = second.representative_hit
    if first_hit is None or second_hit is None:
        return False
    return (
        first_hit.query_start == second_hit.target_start
        and first_hit.query_end == second_hit.target_end
        and first_hit.target_start == second_hit.query_start
        and first_hit.target_end == second_hit.query_end
        and math.isclose(left.identity, right.identity, rel_tol=0.0, abs_tol=1e-12)
    )


def _decision_reasons(decisions: list[PairRelationship]) -> set[str]:
    reasons: set[str] = set()
    for decision in decisions:
        reasons.update(decision.relationship.reasons)
        reasons.update(decision.ambiguity_reasons)
    return reasons


def _edge_from_decision(
    pair: tuple[str, str],
    decision: PairRelationship,
    kind: GraphEdgeKind,
    relationship_type: RelationshipType,
    reasons: tuple[str, ...],
    *,
    retain_coordinates: bool,
    accepted_hit_count: int,
    rejected_hit_count: int,
) -> GraphEdge:
    relationship = decision.relationship
    hit = decision.representative_hit if retain_coordinates else None
    digest = hashlib.sha256(f"{pair[0]}\0{pair[1]}\0{kind.value}".encode()).hexdigest()[:20]
    return GraphEdge(
        edge_id=f"edge_{digest}",
        kind=kind,
        relationship_type=relationship_type,
        query_id=relationship.query_id,
        target_id=relationship.target_id,
        orientation=relationship.orientation,
        query_start=None if hit is None else hit.query_start,
        query_end=None if hit is None else hit.query_end,
        target_start=None if hit is None else hit.target_start,
        target_end=None if hit is None else hit.target_end,
        identity=relationship.identity,
        aligned_length=relationship.aligned_length,
        accepted_hit_count=accepted_hit_count,
        rejected_hit_count=rejected_hit_count,
        reasons=reasons,
    )


def _components(
    node_ids: tuple[str, ...], edges: tuple[GraphEdge, ...]
) -> tuple[MergeComponent, ...]:
    adjacency: dict[str, set[str]] = {identifier: set() for identifier in node_ids}
    edges_by_node: dict[str, list[GraphEdge]] = defaultdict(list)
    for edge in edges:
        adjacency[edge.query_id].add(edge.target_id)
        adjacency[edge.target_id].add(edge.query_id)
        edges_by_node[edge.query_id].append(edge)
        edges_by_node[edge.target_id].append(edge)

    components: list[MergeComponent] = []
    unseen = set(node_ids)
    while unseen:
        start = min(unseen)
        queue = deque((start,))
        members: set[str] = set()
        while queue:
            node = queue.popleft()
            if node in members:
                continue
            members.add(node)
            unseen.discard(node)
            queue.extend(sorted(adjacency[node] - members))
        ordered_members = tuple(sorted(members))
        component_edges = tuple(
            sorted(
                {edge for member in ordered_members for edge in edges_by_node.get(member, ())},
                key=lambda edge: edge.edge_id,
            )
        )
        reasons = _component_ambiguity(ordered_members, component_edges)
        digest = hashlib.sha256("\0".join(ordered_members).encode("utf-8")).hexdigest()[:20]
        components.append(
            MergeComponent(
                component_id=f"component_{digest}",
                sequence_ids=ordered_members,
                relationship_ids=tuple(edge.edge_id for edge in component_edges),
                ambiguous=bool(reasons),
                ambiguity_reasons=reasons,
            )
        )
    return tuple(sorted(components, key=lambda item: item.sequence_ids))


def _component_ambiguity(members: tuple[str, ...], edges: tuple[GraphEdge, ...]) -> tuple[str, ...]:
    reasons: set[str] = set()
    if any(edge.kind is GraphEdgeKind.AMBIGUOUS for edge in edges):
        reasons.add("component contains ambiguous pairwise evidence")

    containment_parents: dict[str, set[str]] = defaultdict(set)
    ports: dict[tuple[str, str], set[str]] = defaultdict(set)
    overlap_edges = [edge for edge in edges if edge.kind is GraphEdgeKind.OVERLAP]
    for edge in edges:
        if edge.kind is GraphEdgeKind.CONTAINMENT:
            contained, container = _contained_and_container(edge)
            containment_parents[contained].add(container)
        elif edge.kind is GraphEdgeKind.OVERLAP:
            for port in _overlap_ports(edge):
                ports[port].add(edge.edge_id)
    if any(len(parents) > 1 for parents in containment_parents.values()):
        reasons.add("contained sequence has multiple possible containers")
    if any(len(edge_ids) > 1 for edge_ids in ports.values()):
        reasons.add("multiple overlaps compete for the same oriented terminal")
    if _has_overlap_cycle(overlap_edges):
        reasons.add("overlap subgraph contains a cycle")
    if _has_orientation_conflict(members, overlap_edges):
        reasons.add("overlap orientations are mutually inconsistent")
    return tuple(sorted(reasons))


def _contained_and_container(edge: GraphEdge) -> tuple[str, str]:
    if edge.relationship_type is RelationshipType.QUERY_CONTAINED_IN_TARGET:
        return edge.query_id, edge.target_id
    return edge.target_id, edge.query_id


def _overlap_ports(edge: GraphEdge) -> tuple[tuple[str, str], tuple[str, str]]:
    target_prefix = "prefix" if edge.orientation is Orientation.FORWARD else "suffix"
    target_suffix = "suffix" if edge.orientation is Orientation.FORWARD else "prefix"
    if edge.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX:
        return (edge.query_id, "suffix"), (edge.target_id, target_prefix)
    return (edge.query_id, "prefix"), (edge.target_id, target_suffix)


def _has_orientation_conflict(members: tuple[str, ...], edges: list[GraphEdge]) -> bool:
    adjacency: dict[str, list[tuple[str, bool]]] = defaultdict(list)
    for edge in edges:
        reverse = edge.orientation is Orientation.REVERSE
        adjacency[edge.query_id].append((edge.target_id, reverse))
        adjacency[edge.target_id].append((edge.query_id, reverse))
    assigned: dict[str, bool] = {}
    for start in members:
        if start in assigned or start not in adjacency:
            continue
        assigned[start] = False
        queue = deque((start,))
        while queue:
            node = queue.popleft()
            for neighbour, reverse in adjacency[node]:
                expected = assigned[node] ^ reverse
                if neighbour in assigned and assigned[neighbour] != expected:
                    return True
                if neighbour not in assigned:
                    assigned[neighbour] = expected
                    queue.append(neighbour)
    return False


def _has_overlap_cycle(edges: list[GraphEdge]) -> bool:
    adjacency: dict[str, list[str]] = defaultdict(list)
    for edge in edges:
        adjacency[edge.query_id].append(edge.target_id)
        adjacency[edge.target_id].append(edge.query_id)
    visited: set[str] = set()
    for start in sorted(adjacency):
        if start in visited:
            continue
        visited.add(start)
        queue = deque(((start, ""),))
        while queue:
            node, parent = queue.popleft()
            for neighbour in adjacency[node]:
                if neighbour == parent:
                    continue
                if neighbour in visited:
                    return True
                visited.add(neighbour)
                queue.append((neighbour, node))
    return False
