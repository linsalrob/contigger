"""Deterministic evaluation of pair classification against checked-in truth."""

from __future__ import annotations

import csv
import hashlib
import json
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, TextIO

from contigger import __version__
from contigger.aligners.minimap2 import parse_paf
from contigger.exceptions import InputValidationError
from contigger.models import (
    AlignmentHit,
    Orientation,
    PairRelationship,
    RelationshipType,
    RunConfig,
)
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.textio import open_text

MERGE_LIKE = frozenset(
    {
        RelationshipType.EXACT_MATCH,
        RelationshipType.QUERY_CONTAINED_IN_TARGET,
        RelationshipType.TARGET_CONTAINED_IN_QUERY,
        RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        RelationshipType.TARGET_SUFFIX_TO_QUERY_PREFIX,
    }
)
TRUTH_COLUMNS = (
    "query_id",
    "target_id",
    "expected_relationship",
    "expected_orientation",
    "expected_status",
    "merge_allowed",
    "expected_overlap_length",
    "expected_identity",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "ambiguity_group",
    "case_id",
    "reason",
)


@dataclass(frozen=True, slots=True)
class ExpectedRelationship:
    """One construction-derived ordered-pair truth record."""

    query_id: str
    target_id: str
    relationship_type: RelationshipType
    orientation: Orientation | None
    status: str
    merge_allowed: bool
    overlap_length: int
    identity: float
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    ambiguity_group: str
    case_id: str
    reason: str


@dataclass(frozen=True, slots=True)
class ObservedRelationship:
    """Pair-classifier output retained independently from benchmark truth."""

    query_id: str
    target_id: str
    relationship_type: RelationshipType
    orientation: Orientation
    status: str
    representative_alignment: AlignmentHit | None
    accepted_hit_count: int
    rejected_hit_count: int
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BenchmarkCaseResult:
    """Scored expected or unexpected ordered pair with explicit diagnostics."""

    query_id: str
    target_id: str
    expected_type: str | None
    expected_orientation: str | None
    observed_type: str | None
    observed_orientation: str | None
    correct: bool
    false_merge: bool
    missed_relationship: bool
    wrong_relationship_type: bool
    wrong_orientation: bool
    ambiguous: bool
    absent_from_paf: bool
    unexpected: bool
    graph_level_deferred: bool
    threshold_boundary_failure: bool
    case_id: str
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BenchmarkSummary:
    """Deterministic aggregate counts and timings for one PAF evaluation."""

    dataset_version: str
    paf_file: str
    expected_relationships: int
    observed_ordered_pairs: int
    self_alignments_excluded: int
    unexpected_observed_pairs: int
    expected_pairs_absent_from_paf: int
    correct_classifications: int
    correct_by_relationship_type: dict[str, int]
    false_merges: int
    missed_relationships: int
    incorrect_relationship_types: int
    incorrect_orientations: int
    pair_level_ambiguous_results: int
    graph_level_ambiguity_cases_deferred: int
    threshold_boundary_failures: int
    observed_by_relationship_type: dict[str, int]
    elapsed_parsing_seconds: float
    elapsed_classification_seconds: float
    elapsed_scoring_seconds: float


@dataclass(frozen=True, slots=True)
class BenchmarkReport:
    """Complete deterministic benchmark result plus reproducibility metadata."""

    summary: BenchmarkSummary
    cases: tuple[BenchmarkCaseResult, ...]
    configuration: dict[str, Any]
    metadata_checksum: str
    paf_checksum: str
    minimap2_preset: str
    minimap2_version: str
    contigger_version: str = __version__


def load_truth(path: Path) -> tuple[ExpectedRelationship, ...]:
    """Parse and validate a line-aware benchmark truth TSV."""
    try:
        with open_text(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(TRUTH_COLUMNS):
                raise InputValidationError(f"{path}:1: unexpected truth-table columns")
            rows: list[ExpectedRelationship] = []
            seen: set[tuple[str, str]] = set()
            for line_number, row in enumerate(reader, start=2):
                try:
                    parsed = _parse_truth_row(row)
                except (KeyError, TypeError, ValueError) as error:
                    raise InputValidationError(f"{path}:{line_number}: {error}") from error
                pair = (parsed.query_id, parsed.target_id)
                if pair in seen:
                    raise InputValidationError(
                        f"{path}:{line_number}: duplicate ordered truth pair "
                        f"{pair[0]!r}, {pair[1]!r}"
                    )
                seen.add(pair)
                rows.append(parsed)
    except InputValidationError:
        raise
    return tuple(sorted(rows, key=lambda item: (item.query_id, item.target_id)))


def _parse_truth_row(row: dict[str, str]) -> ExpectedRelationship:
    if not row["query_id"] or not row["target_id"] or not row["case_id"]:
        raise ValueError("query, target, and case identifiers cannot be empty")
    relationship = RelationshipType(row["expected_relationship"])
    orientation_text = row["expected_orientation"]
    orientation = None if orientation_text == "." else Orientation(orientation_text)
    status = row["expected_status"]
    if status not in {"VALID", "FORBIDDEN", "AMBIGUOUS"}:
        raise ValueError(f"invalid expected status {status!r}")
    boolean = row["merge_allowed"]
    if boolean not in {"true", "false"}:
        raise ValueError(f"invalid Boolean {boolean!r}")
    numeric = [int(row[name]) for name in TRUTH_COLUMNS[6:7] + TRUTH_COLUMNS[8:12]]
    if any(value < 0 for value in numeric):
        raise ValueError("truth coordinates and lengths cannot be negative")
    if numeric[1] > numeric[2] or numeric[3] > numeric[4]:
        raise ValueError("truth interval starts cannot exceed interval ends")
    identity = float(row["expected_identity"])
    if not 0 <= identity <= 100:
        raise ValueError("expected identity must be between 0 and 100")
    return ExpectedRelationship(
        query_id=row["query_id"],
        target_id=row["target_id"],
        relationship_type=relationship,
        orientation=orientation,
        status=status,
        merge_allowed=boolean == "true",
        overlap_length=numeric[0],
        identity=identity,
        query_start=numeric[1],
        query_end=numeric[2],
        target_start=numeric[3],
        target_end=numeric[4],
        ambiguity_group=row["ambiguity_group"],
        case_id=row["case_id"],
        reason=row["reason"],
    )


def evaluate_benchmark(dataset: Path, paf: Path, config: RunConfig) -> BenchmarkReport:
    """Validate dataset structure, classify a complete PAF, and score sparse truth."""
    metadata_path = dataset / "metadata" / "benchmark.json"
    truth_path = dataset / "expected" / "expected_relationships.tsv"
    version_path = dataset / "VERSION"
    tools_path = dataset / "metadata" / "tool_versions.tsv"
    for required in (metadata_path, truth_path, version_path, tools_path, paf):
        if not required.is_file():
            raise InputValidationError(f"benchmark dataset is missing required file: {required}")
    with open_text(metadata_path) as handle:
        try:
            metadata = json.load(handle)
        except (json.JSONDecodeError, TypeError) as error:
            raise InputValidationError(
                f"invalid benchmark metadata {metadata_path}: {error}"
            ) from error
    with open_text(version_path) as handle:
        version = handle.read().strip()
    if not version or metadata.get("benchmark_version") != version:
        raise InputValidationError("benchmark VERSION does not match metadata benchmark_version")
    truth = load_truth(truth_path)
    if metadata.get("relationship_counts") and sum(metadata["relationship_counts"].values()) != len(
        truth
    ):
        raise InputValidationError(
            "benchmark metadata relationship count does not match truth table"
        )

    parsing_start = time.perf_counter()
    self_count = 0
    hits: list[AlignmentHit] = []
    with open_text(paf) as handle:
        for hit in parse_paf(handle):
            if hit.query_id == hit.target_id:
                self_count += 1
            else:
                hits.append(hit)
    parsing_elapsed = time.perf_counter() - parsing_start
    classification_start = time.perf_counter()
    decisions = tuple(classify_pair(group, config) for group in group_ordered_pairs(hits))
    classification_elapsed = time.perf_counter() - classification_start
    scoring_start = time.perf_counter()
    report = _score(paf, config, metadata_path, tools_path, version, truth, decisions, self_count)
    scoring_elapsed = time.perf_counter() - scoring_start
    summary = report.summary
    return BenchmarkReport(
        BenchmarkSummary(
            **{
                **asdict(summary),
                "elapsed_parsing_seconds": parsing_elapsed,
                "elapsed_classification_seconds": classification_elapsed,
                "elapsed_scoring_seconds": scoring_elapsed,
            }
        ),
        report.cases,
        report.configuration,
        report.metadata_checksum,
        report.paf_checksum,
        report.minimap2_preset,
        report.minimap2_version,
    )


def _observed(decision: PairRelationship) -> ObservedRelationship:
    rel = decision.relationship
    return ObservedRelationship(
        rel.query_id,
        rel.target_id,
        rel.relationship_type,
        rel.orientation,
        rel.status,
        decision.representative_hit,
        len(decision.accepted_hits),
        len(decision.rejected_alignments),
        tuple(sorted(set(rel.reasons) | set(decision.ambiguity_reasons))),
    )


def _score(
    paf: Path,
    config: RunConfig,
    metadata_path: Path,
    tools_path: Path,
    version: str,
    truth: tuple[ExpectedRelationship, ...],
    decisions: tuple[PairRelationship, ...],
    self_count: int,
) -> BenchmarkReport:
    expected = {(item.query_id, item.target_id): item for item in truth}
    observed = {(item.query_id, item.target_id): item for item in map(_observed, decisions)}
    graph_groups = {item.ambiguity_group for item in truth if item.ambiguity_group}
    cases = [_score_expected(item, observed.get((item.query_id, item.target_id))) for item in truth]
    for pair in sorted(set(observed) - set(expected)):
        item = observed[pair]
        cases.append(
            BenchmarkCaseResult(
                *pair,
                None,
                None,
                item.relationship_type.value,
                item.orientation.value,
                False,
                False,
                False,
                False,
                False,
                item.relationship_type is RelationshipType.AMBIGUOUS_OVERLAP,
                False,
                True,
                False,
                False,
                "",
                ("ordered pair is absent from sparse truth table",),
            )
        )
    cases.sort(key=lambda item: (not item.false_merge, item.query_id, item.target_id))
    truth_cases = [item for item in cases if not item.unexpected]
    correct_counts = Counter(
        item.expected_type
        for item in truth_cases
        if item.correct and item.expected_type is not None
    )
    observed_counts = Counter(item.relationship_type.value for item in observed.values())
    incorrect_types = sum(item.wrong_relationship_type for item in truth_cases)
    incorrect_orientations = sum(item.wrong_orientation for item in truth_cases)
    summary = BenchmarkSummary(
        version,
        str(paf),
        len(truth),
        len(observed),
        self_count,
        sum(item.unexpected for item in cases),
        sum(item.absent_from_paf for item in truth_cases),
        sum(item.correct for item in truth_cases),
        dict(sorted(correct_counts.items())),
        sum(item.false_merge for item in truth_cases),
        sum(item.missed_relationship for item in truth_cases),
        incorrect_types,
        incorrect_orientations,
        sum(item.ambiguous and not item.graph_level_deferred for item in truth_cases),
        len(graph_groups),
        sum(item.threshold_boundary_failure for item in truth_cases),
        dict(sorted(observed_counts.items())),
        0.0,
        0.0,
        0.0,
    )
    versions = _tool_versions(tools_path)
    preset = "asm5" if "asm5" in paf.name else "asm20" if "asm20" in paf.name else "unknown"
    return BenchmarkReport(
        summary,
        tuple(cases),
        config.as_dict(),
        _sha256(metadata_path),
        _sha256(paf),
        preset,
        versions.get("minimap2", "unknown"),
    )


def _score_expected(
    expected: ExpectedRelationship, observed: ObservedRelationship | None
) -> BenchmarkCaseResult:
    graph_deferred = bool(expected.ambiguity_group)
    absent = observed is None
    observed_type = observed.relationship_type if observed else None
    merge_observed = observed_type in MERGE_LIKE if observed_type else False
    false_merge = not graph_deferred and not expected.merge_allowed and merge_observed
    unambiguous_valid = expected.status == "VALID" and expected.relationship_type in MERGE_LIKE
    missed = unambiguous_valid and (absent or observed_type is RelationshipType.NO_RELATIONSHIP)
    orientation_wrong = (
        observed is not None
        and expected.orientation is not None
        and observed.orientation is not expected.orientation
    )
    expected_observation = (
        expected.relationship_type if expected.merge_allowed else RelationshipType.NO_RELATIONSHIP
    )
    correct = (
        not graph_deferred
        and observed is not None
        and observed_type is expected_observation
        and (expected_observation is RelationshipType.NO_RELATIONSHIP or not orientation_wrong)
    )
    ambiguous = observed_type is RelationshipType.AMBIGUOUS_OVERLAP
    wrong_type = (
        observed is not None
        and observed_type is not expected_observation
        and not graph_deferred
        and not ambiguous
    )
    wrong_orientation = orientation_wrong and not graph_deferred and not ambiguous
    boundary = expected.case_id.startswith(
        ("identity_", "length_", "end_tolerance_", "containment_coverage_")
    )
    reasons: list[str] = []
    if graph_deferred:
        reasons.append("graph-level ambiguity requires multi-target or component context")
    if absent:
        reasons.append("expected ordered pair is absent from PAF")
    if false_merge:
        reasons.append("merge-like pairwise result is forbidden by truth")
    if missed:
        reasons.append("unambiguous valid relationship was missed")
    if observed and observed_type != expected_observation and not graph_deferred:
        reasons.append("observed relationship type differs from truth")
    if orientation_wrong and not graph_deferred:
        reasons.append("observed orientation differs from truth")
    return BenchmarkCaseResult(
        expected.query_id,
        expected.target_id,
        expected.relationship_type.value,
        None if expected.orientation is None else expected.orientation.value,
        None if observed_type is None else observed_type.value,
        None if observed is None else observed.orientation.value,
        correct,
        false_merge,
        missed,
        wrong_type,
        wrong_orientation,
        ambiguous,
        absent,
        False,
        graph_deferred,
        boundary and not (correct or graph_deferred),
        expected.case_id,
        tuple(reasons),
    )


def write_json(report: BenchmarkReport, output: TextIO) -> None:
    """Write a stable JSON report, omitting volatile elapsed timings."""
    payload = asdict(report)
    for name in tuple(payload["summary"]):
        if name.startswith("elapsed_"):
            payload["summary"].pop(name)
    json.dump(payload, output, indent=2, sort_keys=True)
    output.write("\n")


def write_tsv(report: BenchmarkReport, output: TextIO) -> None:
    """Write deterministic per-pair benchmark details."""
    fields = tuple(BenchmarkCaseResult.__dataclass_fields__)
    output.write("\t".join(fields) + "\n")
    for case in report.cases:
        values = asdict(case)
        output.write(
            "\t".join(
                "; ".join(value)
                if isinstance(value, tuple)
                else str(value).lower()
                if isinstance(value, bool)
                else ""
                if value is None
                else str(value)
                for value in (values[field] for field in fields)
            )
            + "\n"
        )


def format_summary(report: BenchmarkReport) -> str:
    """Format a false-merge-first human-readable report."""
    s = report.summary
    false_details = [
        f"  {c.query_id} -> {c.target_id}: {c.observed_type}" for c in report.cases if c.false_merge
    ]
    lines = [
        f"false merges: {s.false_merges}",
        *false_details,
        f"dataset version: {s.dataset_version}",
        f"PAF file: {s.paf_file}",
        "configuration thresholds: " + json.dumps(report.configuration, sort_keys=True),
        f"expected relationships: {s.expected_relationships}",
        f"observed ordered pairs: {s.observed_ordered_pairs}",
        f"self-alignments excluded: {s.self_alignments_excluded}",
        f"unexpected observed pairs: {s.unexpected_observed_pairs}",
        f"expected pairs absent from PAF: {s.expected_pairs_absent_from_paf}",
        f"correct classifications: {s.correct_classifications}",
        "correct classifications by relationship type: "
        + json.dumps(s.correct_by_relationship_type, sort_keys=True),
        f"missed relationships: {s.missed_relationships}",
        f"incorrect relationship types: {s.incorrect_relationship_types}",
        f"incorrect orientations: {s.incorrect_orientations}",
        f"pair-level ambiguous results: {s.pair_level_ambiguous_results}",
        f"graph-level ambiguity cases deferred: {s.graph_level_ambiguity_cases_deferred}",
        f"threshold-boundary failures: {s.threshold_boundary_failures}",
        f"elapsed parsing time: {s.elapsed_parsing_seconds:.6f} s",
        f"elapsed classification time: {s.elapsed_classification_seconds:.6f} s",
        f"elapsed scoring time: {s.elapsed_scoring_seconds:.6f} s",
    ]
    return "\n".join(lines)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tool_versions(path: Path) -> dict[str, str]:
    with open_text(path, newline="") as handle:
        return {row["tool"]: row["version"] for row in csv.DictReader(handle, delimiter="\t")}
