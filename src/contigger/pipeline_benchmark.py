"""Stage-aware benchmark for catalogue, candidates, and pair classification."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO

from contigger.benchmark import MERGE_LIKE, ExpectedRelationship, evaluate_benchmark, load_truth
from contigger.catalogue import build_catalogue, load_source_sequences
from contigger.exceptions import InputValidationError
from contigger.manifest import parse_manifest
from contigger.minimisers import generate_candidates
from contigger.models import CatalogueMember, Orientation, RelationshipType, RunConfig


@dataclass(frozen=True, slots=True)
class PipelineBenchmarkSummary:
    """Deterministic stage-level recall and safety counts."""

    dataset_version: str
    paf_file: str
    source_contigs: int
    canonical_sequences: int
    possible_canonical_pairs: int
    candidate_pairs: int
    exact_truth_rows: int
    exact_truth_rows_recovered: int
    valid_nonexact_case_groups: int
    valid_nonexact_case_groups_recalled: int
    valid_nonexact_truth_rows: int
    valid_nonexact_truth_rows_recalled: int
    candidate_stage_missed_relationships: int
    relationship_stage_false_merges: int
    relationship_stage_missed_relationships: int
    graph_level_ambiguity_cases_deferred: int


@dataclass(frozen=True, slots=True)
class PipelineBenchmarkReport:
    """A benchmark report joining deterministic results across pipeline stages."""

    summary: PipelineBenchmarkSummary
    configuration: dict[str, object]
    metadata_checksum: str
    paf_checksum: str
    minimap2_preset: str
    minimap2_version: str


def evaluate_pipeline_benchmark(
    dataset: Path,
    paf: Path,
    config: RunConfig,
    *,
    kmer_size: int = 21,
    window_size: int = 10,
    min_shared_minimisers: int = 5,
    max_minimiser_frequency: int = 100,
    terminal_band: int = 1000,
) -> PipelineBenchmarkReport:
    """Evaluate exact recovery, candidate recall, then eligible pair decisions.

    Checked-in complete PAFs provide deterministic alignment observations; this
    command does not invoke minimap2 and does not construct biological merges.
    """
    validation = parse_manifest(dataset / "manifest.tsv")
    catalogue = build_catalogue(load_source_sequences(validation.samples))
    truth = load_truth(dataset / "expected" / "expected_relationships.tsv")
    pair_report = evaluate_benchmark(dataset, paf, config)
    candidates = generate_candidates(
        catalogue.sequences,
        kmer_size=kmer_size,
        window_size=window_size,
        min_shared_minimisers=min_shared_minimisers,
        max_minimiser_frequency=max_minimiser_frequency,
        terminal_band=terminal_band,
    )
    candidate_pairs = {(item.query_id, item.target_id) for item in candidates}
    members = _source_members(catalogue.members)
    case_results = {
        (item.query_id, item.target_id): item for item in pair_report.cases if not item.unexpected
    }

    exact = [item for item in truth if item.relationship_type is RelationshipType.EXACT_MATCH]
    exact_recovered = sum(_exact_recovered(item, members) for item in exact)
    valid_nonexact = [
        item
        for item in truth
        if item.status == "VALID"
        and item.merge_allowed
        and item.relationship_type in MERGE_LIKE
        and item.relationship_type is not RelationshipType.EXACT_MATCH
    ]
    recalled_rows = [
        item for item in valid_nonexact if _catalogue_pair(item, members) in candidate_pairs
    ]
    expected_groups = {item.case_id for item in valid_nonexact}
    recalled_groups = {item.case_id for item in recalled_rows}
    eligible_truth = [
        item
        for item in truth
        if item.relationship_type is not RelationshipType.EXACT_MATCH
        and _catalogue_pair(item, members) in candidate_pairs
    ]
    try:
        eligible_results = [case_results[(item.query_id, item.target_id)] for item in eligible_truth]
    except KeyError as error:
        raise InputValidationError(
            f"benchmark truth references a pair with no evaluation result: {error.args[0]}"
        ) from error
    configuration: dict[str, object] = {
        **config.as_dict(),
        "kmer_size": kmer_size,
        "window_size": window_size,
        "min_shared_minimisers": min_shared_minimisers,
        "max_minimiser_frequency": max_minimiser_frequency,
        "terminal_band": terminal_band,
    }
    canonical_count = len(catalogue.sequences)
    summary = PipelineBenchmarkSummary(
        dataset_version=pair_report.summary.dataset_version,
        paf_file=str(paf),
        source_contigs=len(catalogue.members),
        canonical_sequences=canonical_count,
        possible_canonical_pairs=canonical_count * (canonical_count - 1) // 2,
        candidate_pairs=len(candidates),
        exact_truth_rows=len(exact),
        exact_truth_rows_recovered=exact_recovered,
        valid_nonexact_case_groups=len(expected_groups),
        valid_nonexact_case_groups_recalled=len(recalled_groups),
        valid_nonexact_truth_rows=len(valid_nonexact),
        valid_nonexact_truth_rows_recalled=len(recalled_rows),
        candidate_stage_missed_relationships=len(valid_nonexact) - len(recalled_rows),
        relationship_stage_false_merges=sum(item.false_merge for item in eligible_results),
        relationship_stage_missed_relationships=(
            len(valid_nonexact)
            - len(recalled_rows)
            + sum(
                case_results[(item.query_id, item.target_id)].missed_relationship
                for item in recalled_rows
            )
        ),
        graph_level_ambiguity_cases_deferred=(
            pair_report.summary.graph_level_ambiguity_cases_deferred
        ),
    )
    return PipelineBenchmarkReport(
        summary,
        configuration,
        pair_report.metadata_checksum,
        pair_report.paf_checksum,
        pair_report.minimap2_preset,
        pair_report.minimap2_version,
    )


def format_pipeline_summary(report: PipelineBenchmarkReport) -> str:
    """Format a false-merge-first pipeline summary."""
    summary = report.summary
    return "\n".join(
        (
            f"relationship-stage false merges: {summary.relationship_stage_false_merges}",
            f"dataset version: {summary.dataset_version}",
            f"PAF file: {summary.paf_file}",
            f"source contigs: {summary.source_contigs}",
            f"canonical sequences: {summary.canonical_sequences}",
            f"possible canonical pairs: {summary.possible_canonical_pairs}",
            f"candidate pairs: {summary.candidate_pairs}",
            "exact truth rows recovered: "
            f"{summary.exact_truth_rows_recovered}/{summary.exact_truth_rows}",
            "valid non-exact case groups recalled: "
            f"{summary.valid_nonexact_case_groups_recalled}/{summary.valid_nonexact_case_groups}",
            f"candidate-stage missed relationships: {summary.candidate_stage_missed_relationships}",
            "relationship-stage missed relationships: "
            f"{summary.relationship_stage_missed_relationships}",
            f"graph-level ambiguity cases deferred: {summary.graph_level_ambiguity_cases_deferred}",
            "no graph was constructed and no contigs were merged",
        )
    )


def write_pipeline_json(report: PipelineBenchmarkReport, output: TextIO) -> None:
    """Write deterministic pipeline benchmark JSON."""
    json.dump(asdict(report), output, indent=2, sort_keys=True)
    output.write("\n")


def _source_members(members: tuple[CatalogueMember, ...]) -> dict[str, CatalogueMember]:
    result: dict[str, CatalogueMember] = {}
    for member in members:
        original = member.original_identifier
        if original in result:
            raise InputValidationError(f"benchmark source identifier is not unique: {original}")
        result[original] = member
    return result


def _catalogue_pair(
    expected: ExpectedRelationship, members: dict[str, CatalogueMember]
) -> tuple[str, str]:
    try:
        query = members[expected.query_id]
        target = members[expected.target_id]
    except KeyError as error:
        raise InputValidationError(
            f"benchmark truth references unknown source contig: {error.args[0]}"
        ) from error
    first, second = sorted((query.catalogue_id, target.catalogue_id))
    return first, second


def _exact_recovered(expected: ExpectedRelationship, members: dict[str, CatalogueMember]) -> bool:
    try:
        query = members[expected.query_id]
        target = members[expected.target_id]
    except KeyError as error:
        raise InputValidationError(
            f"benchmark truth references unknown source contig: {error.args[0]}"
        ) from error
    if query.catalogue_id != target.catalogue_id:
        return False
    observed = (
        Orientation.FORWARD if query.orientation is target.orientation else Orientation.REVERSE
    )
    return observed is expected.orientation
