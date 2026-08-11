"""Conservative end-to-end catalogue, graph, and sequence construction."""

from __future__ import annotations

import json
import os
import platform
import time
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType

from contigger.aligners.minimap2 import Minimap2Aligner
from contigger.alignment_planning import (
    execute_indexed_selective_alignments,
    plan_selective_alignments,
)
from contigger.catalogue import build_catalogue, load_source_sequences
from contigger.decision_policy import evaluate_graph_decisions
from contigger.evidence.bam import BamEvidenceProvider
from contigger.exceptions import InputValidationError
from contigger.graph import build_relationship_graph
from contigger.minimisers import CandidateGenerationMetrics, generate_candidates_with_metrics
from contigger.models import (
    AlignmentHit,
    AlignmentRequest,
    CandidatePair,
    CatalogueMember,
    CatalogueSequence,
    ContainmentDecision,
    GraphDecisionPlan,
    GraphDecisionStatus,
    GraphEdge,
    LinearPathPlan,
    Orientation,
    PairRelationship,
    PathPlanningResult,
    RelationshipGraph,
    RelationshipType,
    RunConfig,
    SampleInput,
    SequenceCatalogue,
    SequenceRecord,
)
from contigger.outputs import OutputPaths, output_paths
from contigger.path_planning import plan_linear_paths
from contigger.provenance import ProvenanceRecord, write_provenance
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.utilities.sequences import reverse_complement

_resource: ModuleType | None
try:
    import resource as _resource
except ImportError:  # pragma: no cover - exercised on Windows
    _resource = None


def merge_samples(samples: tuple[SampleInput, ...], config: RunConfig) -> tuple[Path, ...]:
    """Run the conservative production pipeline and write deterministic outputs."""
    started = time.monotonic()
    if config.evidence.value == "reads":
        raise InputValidationError(
            "evidence mode 'reads' is not implemented for merge; use 'none' or 'alignments'"
        )
    if config.evidence.value == "alignments":
        missing = tuple(sample.sample for sample in samples if sample.bam is None)
        if missing:
            raise InputValidationError(
                "evidence mode 'alignments' requires a BAM/CRAM for every sample; "
                f"missing: {', '.join(missing)}"
            )
        samtools_versions: dict[str, str | None] = {}
        samtools_commands: dict[str, tuple[tuple[str, ...], ...]] = {}
        for sample in samples:
            provider = BamEvidenceProvider(sample)
            provider.validate_source()
            samtools_versions[sample.sample] = provider.version
            samtools_commands[sample.sample] = provider.commands
    else:
        samtools_versions = {}
        samtools_commands = {}
    stages: dict[str, float] = {}

    stage_start = time.monotonic()
    records = load_source_sequences(samples)
    catalogue = build_catalogue(records)
    stages["catalogue"] = time.monotonic() - stage_start
    print(
        f"Loaded {len(records)} contigs; canonical catalogue {len(catalogue.sequences)} sequences"
    )

    stage_start = time.monotonic()
    candidates, candidate_metrics = generate_candidates_with_metrics(
        catalogue.sequences,
        kmer_size=config.kmer_size,
        window_size=config.window_size,
        min_shared_minimisers=config.min_shared_minimisers,
        max_minimiser_frequency=config.max_minimiser_frequency,
        terminal_band=max(config.min_overlap, config.min_containment),
    )
    if config.max_candidate_pairs is not None and len(candidates) > config.max_candidate_pairs:
        raise InputValidationError(
            f"candidate pair count {len(candidates)} exceeds "
            f"--max-candidate-pairs {config.max_candidate_pairs}; "
            "increase the limit or tighten candidate-generation parameters"
        )
    requests = plan_selective_alignments(catalogue.sequences, candidates)
    stages["candidates"] = time.monotonic() - stage_start
    print(
        f"Generated {len(candidates)} candidate pairs from "
        f"{candidate_metrics.retained_observations} retained minimiser observations"
    )

    stage_start = time.monotonic()
    hits: tuple[AlignmentHit, ...]
    tool_versions: dict[str, str | None] = {"minimap2": None}
    aligner_metrics: dict[str, int] = {"index_builds": 0, "index_reuses": 0, "alignment_batches": 0}
    if requests:
        index_dir = config.index_dir or config.output_prefix.parent / (
            f".{config.output_prefix.name}-indexes"
        )
        aligner = Minimap2Aligner(threads=config.threads, preset=config.minimap2_preset)
        hits = execute_indexed_selective_alignments(requests, aligner, index_dir)
        tool_versions = {"minimap2": aligner.tool_version}
        for name in aligner_metrics:
            aligner_metrics[name] = int(getattr(aligner, name, 0))
    else:
        hits = ()
    stages["alignment"] = time.monotonic() - stage_start
    print(f"Classified alignment input from {len(hits)} alignment observations")

    stage_start = time.monotonic()
    relationships = tuple(classify_pair(group, config) for group in group_ordered_pairs(hits))
    graph = build_relationship_graph(
        relationships, sequence_ids=(item.identifier for item in catalogue.sequences)
    )
    safe_edges = tuple(
        edge.edge_id
        for edge in graph.overlap_edges
        if _edge_has_exact_reconcilable_overlap(edge, catalogue, config.end_tolerance)
    )
    decisions = evaluate_graph_decisions(graph, intrinsically_safe_edge_ids=safe_edges)
    paths = plan_linear_paths(
        catalogue,
        graph,
        intrinsically_safe_edge_ids=safe_edges,
    )
    stages["graph"] = time.monotonic() - stage_start

    stage_start = time.monotonic()
    constructed, provenance, ambiguous, merge_stats = _construct_outputs(
        catalogue, graph, decisions, paths, relationships, config
    )
    stages["construction"] = time.monotonic() - stage_start
    stages["total"] = time.monotonic() - started
    print(
        f"Merged {merge_stats['merged_linear_paths']} safe paths; "
        f"deferred {merge_stats['deferred_junctions']} junctions"
    )

    paths_out = output_paths(config.output_prefix)
    stats = _stats(
        samples,
        records,
        catalogue,
        candidates,
        requests,
        hits,
        relationships,
        graph,
        decisions,
        paths,
        constructed,
        merge_stats,
        config,
        tool_versions,
        stages,
        aligner_metrics,
        samtools_versions,
        samtools_commands,
        candidate_metrics,
    )
    _write_outputs_atomic(
        paths_out,
        constructed,
        provenance,
        relationships,
        ambiguous,
        graph,
        paths,
        stats,
        _join_support_rows(graph, catalogue, config),
        emit_gfa=config.emit_gfa,
    )
    return (
        paths_out.fasta,
        paths_out.provenance,
        paths_out.relationships,
        paths_out.ambiguous,
        paths_out.gfa,
        paths_out.stats,
    )


def _edge_has_exact_reconcilable_overlap(
    edge: GraphEdge, catalogue: SequenceCatalogue, end_tolerance: int
) -> bool:
    """Return whether an overlap is terminal and nucleotide-identical without reads."""
    if edge.query_start is None or edge.query_end is None:
        return False
    if edge.target_start is None or edge.target_end is None:
        return False
    if edge.relationship_type not in {
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX,
    }:
        return False
    lookup = {item.identifier: item for item in catalogue.sequences}
    query = lookup[edge.query_id]
    target = lookup[edge.target_id]
    query_span = edge.query_end - edge.query_start
    target_span = edge.target_end - edge.target_start
    if query_span <= 0 or query_span != target_span or edge.aligned_length != query_span:
        return False
    query_part = query.sequence[edge.query_start : edge.query_end]
    target_sequence = target.sequence
    if edge.orientation is Orientation.REVERSE:
        target_sequence = reverse_complement(target_sequence)
        target_start = target.length - edge.target_end
        target_end = target.length - edge.target_start
    else:
        target_start = edge.target_start
        target_end = edge.target_end
    target_part = target_sequence[target_start:target_end]
    if edge.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX:
        terminal_distances = (query.length - edge.query_end, target_start)
    else:
        terminal_distances = (edge.query_start, target.length - target_end)
    if any(distance > end_tolerance for distance in terminal_distances):
        return False
    return query_part == target_part


def _construct_outputs(
    catalogue: SequenceCatalogue,
    graph: RelationshipGraph,
    decisions: GraphDecisionPlan,
    paths: PathPlanningResult,
    relationships: tuple[PairRelationship, ...],
    config: RunConfig,
) -> tuple[
    tuple[SequenceRecord, ...],
    tuple[ProvenanceRecord, ...],
    tuple[tuple[str, str], ...],
    dict[str, int],
]:
    """Construct safe paths and retain every deferred sequence with diagnostics."""
    lookup = {item.identifier: item for item in catalogue.sequences}
    containment = {
        item.contained_sequence_id: item
        for item in decisions.containment_decisions
        if item.status is GraphDecisionStatus.ELIGIBLE
    }
    constructed: dict[str, SequenceRecord] = {}
    provenance: list[ProvenanceRecord] = []
    ambiguous: list[tuple[str, str]] = []
    removed: set[str] = set()
    path_nodes: set[str] = set()
    successful_paths = 0
    successful_junctions = 0
    edge_lookup = {
        edge.edge_id: edge
        for edge in (*graph.overlap_edges, *graph.containment_edges, *graph.ambiguous_edges)
    }
    for path in paths.paths:
        try:
            merged_sequence, intervals = _construct_path(path, edge_lookup, lookup)
        except InputValidationError as error:
            ambiguous.append((path.path_id, str(error)))
            continue
        output_id = f"merge_{path.path_id.removeprefix('path_')}"
        constructed[output_id] = SequenceRecord(
            output_id,
            "",
            output_id,
            "conservative conflict-free terminal-overlap merge",
            merged_sequence,
            len(merged_sequence),
        )
        successful_paths += 1
        successful_junctions += len(path.edge_ids)
        path_nodes.update(node.sequence_id for node in path.nodes)
        for node, (start, end) in zip(path.nodes, intervals, strict=True):
            for path_member in node.source_members:
                source_length = lookup[node.sequence_id].length
                provenance.append(
                    ProvenanceRecord(
                        output_id,
                        path_member.source_sample,
                        path_member.original_identifier,
                        "TERMINAL_OVERLAP",
                        path_member.orientation,
                        0,
                        source_length,
                        start,
                        end,
                        1.0,
                        "constructed_path",
                        "exact terminal overlap; no evidence required",
                        config.evidence.value,
                    )
                )
    for catalogue_sequence in catalogue.sequences:
        if catalogue_sequence.identifier in path_nodes:
            continue
        if catalogue_sequence.identifier in containment:
            removed.add(catalogue_sequence.identifier)
            continue
        constructed[catalogue_sequence.identifier] = SequenceRecord(
            catalogue_sequence.identifier,
            "",
            catalogue_sequence.identifier,
            "canonical catalogue representative",
            catalogue_sequence.sequence,
            catalogue_sequence.length,
        )
        for member in _members(catalogue, catalogue_sequence.identifier):
            provenance.append(
                _catalogue_record(
                    member,
                    catalogue_sequence.identifier,
                    catalogue_sequence.length,
                    "canonical representative" if member.representative else "exact duplicate",
                    "exact or reverse-complement catalogue identity",
                    config.evidence.value,
                )
            )
    for contained_id, decision in containment.items():
        container, output_start, output_end = _containment_root(
            contained_id, containment, edge_lookup, lookup
        )
        if container not in constructed:
            ambiguous.append((contained_id, "eligible containment container was not emitted"))
            continue
        for member in _members(catalogue, contained_id):
            provenance.append(
                _catalogue_record(
                    member,
                    container,
                    lookup[contained_id].length,
                    "contained_removed",
                    "; ".join(decision.reasons),
                    config.evidence.value,
                    output_start=output_start,
                    output_end=output_end,
                )
            )
    deferred_ids = {component_id for component_id in paths.deferred_component_ids}
    for component in graph.components:
        if component.component_id in deferred_ids or component.ambiguous:
            ambiguous.append(
                (
                    component.component_id,
                    "; ".join(component.ambiguity_reasons)
                    or "component deferred by conservative policy",
                )
            )
    for relationship in relationships:
        if relationship.relationship.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP:
            ambiguous.append(
                (
                    f"{relationship.relationship.query_id}:{relationship.relationship.target_id}",
                    "; ".join(relationship.relationship.reasons),
                )
            )
    merge_stats = {
        "contained_contigs_removed": len(removed),
        "merged_linear_paths": successful_paths,
        "merged_junctions": successful_junctions,
        "deferred_junctions": len(graph.overlap_edges) - successful_junctions,
    }
    return (
        tuple(sorted(constructed.values(), key=lambda item: item.identifier)),
        tuple(provenance),
        tuple(sorted(set(ambiguous))),
        merge_stats,
    )


def _construct_path(
    path: LinearPathPlan,
    edges: dict[str, GraphEdge],
    sequences: dict[str, CatalogueSequence],
) -> tuple[str, tuple[tuple[int, int], ...]]:
    """Construct one path, refusing any coordinate or base disagreement."""
    first = path.nodes[0]
    sequence = _oriented(sequences[first.sequence_id].sequence, first.orientation)
    intervals: list[tuple[int, int]] = [(0, len(sequence))]
    for index, edge_id in enumerate(path.edge_ids):
        edge = edges[edge_id]
        current = path.nodes[index]
        next_node = path.nodes[index + 1]
        next_sequence = _oriented(sequences[next_node.sequence_id].sequence, next_node.orientation)
        hit_query = sequences[edge.query_id]
        hit_target = sequences[edge.target_id]
        if (
            edge.query_start is None
            or edge.query_end is None
            or edge.target_start is None
            or edge.target_end is None
        ):
            raise InputValidationError(f"edge {edge_id} lacks construction coordinates")
        query_part = _interval(
            hit_query.sequence,
            edge.query_start,
            edge.query_end,
            current.orientation if current.sequence_id == edge.query_id else next_node.orientation,
        )
        target_part = _interval(
            hit_target.sequence,
            edge.target_start,
            edge.target_end,
            current.orientation if current.sequence_id == edge.target_id else next_node.orientation,
        )
        if current.sequence_id == edge.query_id:
            current_part, next_part = query_part, target_part
        else:
            current_part, next_part = target_part, query_part
        overlap = len(current_part)
        if overlap == 0 or len(current_part) != len(next_part) or current_part != next_part:
            raise InputValidationError(
                f"edge {edge_id} has an imperfect or orientation-inconsistent overlap"
            )
        current_start = (
            edge.query_start if current.sequence_id == edge.query_id else edge.target_start
        )
        if current_start is None:
            raise InputValidationError(f"edge {edge_id} lacks current coordinates")
        if (
            not sequence.endswith(current_part)
            or _interval(
                sequences[current.sequence_id].sequence,
                current_start,
                current_start + len(current_part),
                current.orientation,
            )
            != current_part
            or not next_sequence.startswith(next_part)
        ):
            raise InputValidationError(f"edge {edge_id} is not terminal in the planned orientation")
        start = len(sequence) - overlap
        sequence += next_sequence[overlap:]
        intervals.append((start, start + len(next_sequence)))
    return sequence, tuple(intervals)


def _interval(sequence: str, start: int, end: int, orientation: Orientation) -> str:
    """Return a canonical interval expressed in an explicitly oriented string."""
    if orientation is Orientation.FORWARD:
        return sequence[start:end]
    return reverse_complement(sequence[start:end])


def _oriented(sequence: str, orientation: Orientation) -> str:
    return sequence if orientation is Orientation.FORWARD else reverse_complement(sequence)


def _members(catalogue: SequenceCatalogue, identifier: str) -> tuple[CatalogueMember, ...]:
    return tuple(item for item in catalogue.members if item.catalogue_id == identifier)


def _containment_root(
    contained_id: str,
    containment: dict[str, ContainmentDecision],
    edges: dict[str, GraphEdge],
    sequences: dict[str, CatalogueSequence],
) -> tuple[str, int, int]:
    """Resolve an eligible containment chain to its emitted ancestor coordinates."""
    current = contained_id
    start = 0
    end = sequences[contained_id].length
    visited: set[str] = set()
    while current in containment:
        if current in visited:
            raise InputValidationError(f"containment cycle prevents disposition: {current}")
        visited.add(current)
        decision = containment[current]
        edge = edges[decision.edge_id]
        child, parent, child_start, child_end, parent_start, parent_end = _containment_intervals(
            edge
        )
        if child != current:
            raise InputValidationError(f"containment decision does not match child {current}")
        parent_length = sequences[parent].length
        parent_span = parent_end - parent_start
        start = parent_start + (start * parent_span) // parent_length
        end = parent_start + (end * parent_span) // parent_length
        current = parent
    return current, start, end


def _containment_intervals(edge: GraphEdge) -> tuple[str, str, int, int, int, int]:
    """Return child/parent identifiers and aligned intervals for a containment edge."""
    if edge.query_start is None or edge.query_end is None:
        raise InputValidationError(f"containment edge lacks query coordinates: {edge.edge_id}")
    if edge.target_start is None or edge.target_end is None:
        raise InputValidationError(f"containment edge lacks target coordinates: {edge.edge_id}")
    if edge.relationship_type is RelationshipType.QUERY_CONTAINED_IN_TARGET:
        return (
            edge.query_id,
            edge.target_id,
            edge.query_start,
            edge.query_end,
            edge.target_start,
            edge.target_end,
        )
    if edge.relationship_type is RelationshipType.TARGET_CONTAINED_IN_QUERY:
        return (
            edge.target_id,
            edge.query_id,
            edge.target_start,
            edge.target_end,
            edge.query_start,
            edge.query_end,
        )
    raise InputValidationError(f"edge is not a containment relationship: {edge.edge_id}")


def _catalogue_record(
    member: CatalogueMember,
    output_id: str,
    length: int,
    disposition: str,
    reason: str,
    evidence_mode: str,
    *,
    output_start: int = 0,
    output_end: int | None = None,
) -> ProvenanceRecord:
    if output_end is None:
        output_end = length
    return ProvenanceRecord(
        output_id,
        member.source_sample,
        member.original_identifier,
        "EXACT_MATCH" if disposition != "contained_removed" else "CONTAINMENT",
        member.orientation,
        0,
        length,
        output_start,
        output_end,
        1.0,
        disposition,
        reason,
        evidence_mode,
    )


def _stats(
    samples: tuple[SampleInput, ...],
    records: tuple[SequenceRecord, ...],
    catalogue: SequenceCatalogue,
    candidates: tuple[CandidatePair, ...],
    requests: tuple[AlignmentRequest, ...],
    hits: tuple[AlignmentHit, ...],
    relationships: tuple[PairRelationship, ...],
    graph: RelationshipGraph,
    decisions: GraphDecisionPlan,
    paths: PathPlanningResult,
    constructed: tuple[SequenceRecord, ...],
    merge_stats: dict[str, int],
    config: RunConfig,
    tool_versions: dict[str, str | None],
    stages: dict[str, float],
    aligner_metrics: dict[str, int],
    samtools_versions: dict[str, str | None],
    samtools_commands: dict[str, tuple[tuple[str, ...], ...]],
    candidate_metrics: CandidateGenerationMetrics,
) -> dict[str, object]:
    relationship_counts = Counter(
        item.relationship.relationship_type.value for item in relationships
    )
    ordered_samples = tuple(sorted(samples, key=lambda item: item.sample))
    sample_metrics = _sample_metrics(ordered_samples, records)
    representative_orientations = {
        member.catalogue_id: member.orientation
        for member in catalogue.members
        if member.representative
    }
    return {
        "run_status": "completed",
        "input_samples": len(samples),
        "input_contigs": len(records),
        "input_bases": sum(item.length for item in records),
        "input_manifest": [
            {
                "sample": sample.sample,
                "contigs": str(sample.contigs),
                "bam": str(sample.bam) if sample.bam else None,
                "technology": sample.technology,
                "assembly_graph": str(sample.assembly_graph) if sample.assembly_graph else None,
                "metadata": dict(sorted(sample.metadata.items())),
            }
            for sample in ordered_samples
        ],
        "input_by_sample": sample_metrics,
        "canonical_sequences": len(catalogue.sequences),
        "exact_duplicates_collapsed": len(records) - len(catalogue.sequences),
        "reverse_complement_duplicates_collapsed": sum(
            1
            for member in catalogue.members
            if not member.representative
            and member.orientation is not representative_orientations[member.catalogue_id]
        ),
        "reverse_oriented_catalogue_members": sum(
            1 for member in catalogue.members if member.orientation is Orientation.REVERSE
        ),
        "contained_contigs_removed": merge_stats["contained_contigs_removed"],
        "candidate_pairs": len(candidates),
        "candidate_generation": candidate_metrics.as_dict(),
        "alignment_pairs": len(requests),
        "relationship_counts": dict(sorted(relationship_counts.items())),
        "graph_components": len(graph.components),
        "ambiguous_components": sum(item.ambiguous for item in graph.components),
        "eligible_containments": sum(
            item.status is GraphDecisionStatus.ELIGIBLE for item in decisions.containment_decisions
        ),
        "eligible_overlap_components": sum(
            item.status is GraphDecisionStatus.ELIGIBLE for item in decisions.overlap_decisions
        ),
        "merged_linear_paths": merge_stats["merged_linear_paths"],
        "merged_junctions": merge_stats["merged_junctions"],
        "deferred_junctions": merge_stats["deferred_junctions"],
        "output_contigs": len(constructed),
        "output_bases": sum(item.length for item in constructed),
        "tool_versions": tool_versions,
        "samtools_versions": samtools_versions,
        "tool_commands": {
            "samtools": {
                sample: [list(command) for command in commands]
                for sample, commands in sorted(samtools_commands.items())
            }
        },
        "alignment_metrics": aligner_metrics,
        "minimap2_preset": config.minimap2_preset,
        "configuration": config.as_dict(),
        "elapsed_times_by_stage": stages,
        "resource_usage": _resource_usage(),
    }


def _sample_metrics(
    samples: tuple[SampleInput, ...], records: tuple[SequenceRecord, ...]
) -> list[dict[str, object]]:
    """Summarise source contig counts and bases in deterministic sample order."""
    totals: dict[str, list[int]] = {sample.sample: [0, 0] for sample in samples}
    for record in records:
        count_and_bases = totals.setdefault(record.source_sample, [0, 0])
        count_and_bases[0] += 1
        count_and_bases[1] += record.length
    return [
        {
            "sample": sample,
            "contigs": values[0],
            "bases": values[1],
        }
        for sample, values in sorted(totals.items())
    ]


def _resource_usage() -> dict[str, object]:
    """Capture portable process and available Slurm resource metadata."""
    if _resource is None:
        return {"available": False}
    usage = _resource.getrusage(_resource.RUSAGE_SELF)
    # Linux reports KiB; macOS reports bytes. Contigger's production jobs run on
    # Linux, but keeping the conversion explicit makes local baselines comparable.
    peak_rss_kib = int(usage.ru_maxrss)
    if platform.system() == "Darwin":
        peak_rss_kib //= 1024
    resource_usage: dict[str, object] = {
        "available": True,
        "peak_rss_kib": peak_rss_kib,
        "user_cpu_seconds": round(usage.ru_utime, 6),
        "system_cpu_seconds": round(usage.ru_stime, 6),
    }
    slurm_fields = {
        "job_id": "SLURM_JOB_ID",
        "job_partition": "SLURM_JOB_PARTITION",
        "job_nodelist": "SLURM_JOB_NODELIST",
        "cpus_per_task": "SLURM_CPUS_PER_TASK",
        "requested_memory_per_node": "SLURM_MEM_PER_NODE",
        "requested_memory_per_cpu": "SLURM_MEM_PER_CPU",
        "requested_time": "SLURM_TIMELIMIT",
    }
    slurm = {
        name: os.environ[variable]
        for name, variable in slurm_fields.items()
        if os.environ.get(variable)
    }
    if slurm:
        resource_usage["slurm"] = slurm
    return resource_usage


def _write_outputs_atomic(
    paths: OutputPaths,
    sequences: tuple[SequenceRecord, ...],
    provenance: tuple[ProvenanceRecord, ...],
    relationships: tuple[PairRelationship, ...],
    ambiguous: tuple[tuple[str, str], ...],
    graph: RelationshipGraph,
    paths_result: PathPlanningResult,
    stats: dict[str, object],
    join_support_rows: tuple[tuple[str, ...], ...],
    *,
    emit_gfa: bool,
) -> None:
    """Write all requested artifacts through a same-directory temporary set."""
    paths.fasta.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".contigger-output-", dir=paths.fasta.parent) as directory:
        temporary = Path(directory)
        _write_fasta(sequences, temporary / paths.fasta.name)
        write_provenance(temporary / paths.provenance.name, provenance)
        _write_relationships(temporary / paths.relationships.name, relationships)
        _write_ambiguous(temporary / paths.ambiguous.name, ambiguous)
        if emit_gfa:
            _write_gfa(temporary / paths.gfa.name, sequences, graph, paths_result)
        else:
            (temporary / paths.gfa.name).write_text("", encoding="utf-8")
        (temporary / paths.stats.name).write_text(
            json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        _write_join_support(temporary / paths.join_support.name, join_support_rows)
        _write_variants(temporary / paths.variants.name)
        for output in (
            paths.fasta,
            paths.provenance,
            paths.relationships,
            paths.ambiguous,
            paths.gfa,
            paths.join_support,
            paths.variants,
        ):
            (temporary / output.name).replace(output)
        # Install stats last so a failed later replacement can never publish a
        # successful run status for an incomplete output set.
        (temporary / paths.stats.name).replace(paths.stats)


def _write_fasta(sequences: Iterable[SequenceRecord], path: Path) -> None:
    with path.open("w", encoding="ascii", newline="") as handle:
        for record in sequences:
            handle.write(f">{record.identifier}\n")
            for start in range(0, record.length, 80):
                handle.write(record.sequence[start : start + 80] + "\n")


def _write_relationships(path: Path, relationships: Iterable[PairRelationship]) -> None:
    columns = (
        "query",
        "target",
        "relationship_type",
        "orientation",
        "identity",
        "aligned_length",
        "status",
        "reasons",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for item in relationships:
            rel = item.relationship
            handle.write(
                "\t".join(
                    (
                        rel.query_id,
                        rel.target_id,
                        rel.relationship_type.value,
                        rel.orientation.value,
                        f"{rel.identity:.6f}",
                        str(rel.aligned_length),
                        rel.status,
                        "; ".join(rel.reasons),
                    )
                )
                + "\n"
            )


def _write_ambiguous(path: Path, rows: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("component_or_pair\treason\n")
        for identifier, reason in rows:
            handle.write(f"{identifier}\t{reason}\n")


def _join_support_rows(
    graph: RelationshipGraph, catalogue: SequenceCatalogue, config: RunConfig
) -> tuple[tuple[str, ...], ...]:
    """Report imperfect junctions without authorising an unreviewed consensus."""
    if config.evidence.value != "alignments":
        return ()
    rows: list[tuple[str, ...]] = []
    for edge in graph.overlap_edges:
        if _edge_has_exact_reconcilable_overlap(edge, catalogue, config.end_tolerance):
            continue
        rows.append(
            (
                "unknown",
                edge.edge_id,
                edge.query_id,
                edge.target_id,
                "",
                "0",
                "0",
                "0",
                "0",
                "0.0",
                "",
                "",
                "",
                "DEFERRED",
                "",
                "no approved junction policy; imperfect overlap retained",
            )
        )
    return tuple(rows)


def _write_join_support(path: Path, rows: Iterable[tuple[str, ...]]) -> None:
    columns = (
        "component_id",
        "edge_id",
        "left_contig",
        "right_contig",
        "technology",
        "samples_tested",
        "selected_reads",
        "remapped_reads",
        "spanning_reads",
        "spanning_fraction",
        "minimum_mapping_quality",
        "minimum_flank",
        "junction_status",
        "policy_id",
        "decision_reason",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def _write_variants(path: Path) -> None:
    """Write an explicit empty variant report until a reviewed policy resolves sites."""
    path.write_text(
        "edge_id\toverlap_position\tleft_base\tright_base\tdecision\tchosen_base\t"
        "sample\tdepth\tallele_counts\tmean_base_quality\tmean_mapping_quality\t"
        "strand_support\treason\n",
        encoding="utf-8",
    )


def _write_gfa(
    path: Path,
    sequences: Iterable[SequenceRecord],
    graph: RelationshipGraph,
    paths_result: PathPlanningResult,
) -> None:
    ordered_sequences = tuple(sequences)
    emitted = {sequence.identifier for sequence in ordered_sequences}
    with path.open("w", encoding="ascii", newline="") as handle:
        handle.write("H\tVN:Z:1.0\n")
        for sequence in ordered_sequences:
            handle.write(f"S\t{sequence.identifier}\t{sequence.sequence}\n")
        for edge in graph.overlap_edges:
            if edge.query_id not in emitted or edge.target_id not in emitted:
                continue
            target_orientation = "+" if edge.orientation is Orientation.FORWARD else "-"
            handle.write(
                f"L\t{edge.query_id}\t+\t{edge.target_id}\t{target_orientation}\t"
                f"{edge.aligned_length}M\n"
            )
