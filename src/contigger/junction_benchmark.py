"""Typed junction truth, observational scoring, and conservative support policy."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from contigger.exceptions import InputValidationError
from contigger.models import (
    JunctionPolicyReview,
    JunctionSupportDecision,
    JunctionSupportPolicy,
    JunctionSupportStatus,
    TargetedJunctionEvidence,
)
from contigger.textio import open_text

JUNCTION_TRUTH_COLUMNS = (
    "case_id",
    "left_contig",
    "right_contig",
    "junction_is_true",
    "source_coordinate",
    "circular_wrap",
    "selected_spanning_reads",
    "selected_nonspanning_reads",
    "expected_read_evidence",
    "reason",
)

JUNCTION_POLICY_REVIEW_FIELDS = frozenset(
    {
        "truth_dataset_sha256",
        "candidate_baseline_sha256",
        "reviewer",
        "reviewed_at",
        "decision",
        "technology",
        "remapping_preset",
        "minimum_spanning_reads",
        "minimum_spanning_fraction",
        "minimum_spanning_flank",
        "minimum_mapping_quality",
    }
)


@dataclass(frozen=True, slots=True)
class ExpectedJunction:
    """One construction-derived adjacency truth record."""

    case_id: str
    left_contig_id: str
    right_contig_id: str
    junction_is_true: bool
    source_coordinate: int
    circular_wrap: bool
    selected_spanning_reads: int
    selected_nonspanning_reads: int
    expected_read_evidence: str
    reason: str


@dataclass(frozen=True, slots=True)
class JunctionBenchmarkCase:
    """Observed targeted-remapping detection scored against one truth case."""

    case_id: str
    junction_is_true: bool
    spanning_reads: int
    remapped_reads: int
    detected: bool
    correct_detection: bool
    false_support: bool
    missed_support: bool
    absent_observation: bool
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class JunctionBenchmarkSummary:
    """Deterministic aggregate observational counts without merge authorization."""

    expected_cases: int
    true_junctions: int
    artificial_junctions: int
    observed_cases: int
    absent_observations: int
    correct_detections: int
    false_support_cases: int
    missed_support_cases: int


@dataclass(frozen=True, slots=True)
class JunctionBenchmarkReport:
    """Complete ordered junction benchmark result."""

    summary: JunctionBenchmarkSummary
    cases: tuple[JunctionBenchmarkCase, ...]


def load_junction_policy_review(path: Path) -> JunctionPolicyReview:
    """Load one strict JSON review artifact without authorizing a policy."""
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys
        )
    except (OSError, UnicodeError, ValueError) as error:
        raise InputValidationError(f"cannot read junction policy review {path}: {error}") from error
    if not isinstance(payload, dict):
        raise InputValidationError(f"{path}: review artifact must be a JSON object")
    keys = set(payload)
    missing = sorted(JUNCTION_POLICY_REVIEW_FIELDS - keys)
    unknown = sorted(keys - JUNCTION_POLICY_REVIEW_FIELDS - {"notes"})
    if missing:
        raise InputValidationError(f"{path}: review artifact is missing {missing[0]}")
    if unknown:
        raise InputValidationError(f"{path}: review artifact has unknown field {unknown[0]!r}")
    string_fields = {
        "truth_dataset_sha256",
        "candidate_baseline_sha256",
        "reviewer",
        "reviewed_at",
        "decision",
        "technology",
        "remapping_preset",
        "notes",
    }
    for field in string_fields:
        if field in payload and not isinstance(payload[field], str):
            raise InputValidationError(f"{path}: review field {field!r} must be a string")
    for field in ("minimum_spanning_reads", "minimum_spanning_flank", "minimum_mapping_quality"):
        if type(payload[field]) is not int:
            raise InputValidationError(f"{path}: review field {field!r} must be an integer")
    if type(payload["minimum_spanning_fraction"]) not in {int, float}:
        raise InputValidationError(
            f"{path}: review field 'minimum_spanning_fraction' must be a number"
        )
    try:
        return JunctionPolicyReview(
            truth_dataset_sha256=payload["truth_dataset_sha256"],
            candidate_baseline_sha256=payload["candidate_baseline_sha256"],
            reviewer=payload["reviewer"],
            reviewed_at=payload["reviewed_at"],
            decision=payload["decision"],
            technology=payload["technology"],
            remapping_preset=payload["remapping_preset"],
            minimum_spanning_reads=payload["minimum_spanning_reads"],
            minimum_spanning_fraction=payload["minimum_spanning_fraction"],
            minimum_spanning_flank=payload["minimum_spanning_flank"],
            minimum_mapping_quality=payload["minimum_mapping_quality"],
            notes=payload.get("notes", ""),
        )
    except (TypeError, ValueError) as error:
        raise InputValidationError(f"{path}: invalid review artifact value: {error}") from error


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON object members instead of silently choosing one."""
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field {key!r}")
        result[key] = value
    return result


def load_junction_truth(path: Path) -> tuple[ExpectedJunction, ...]:
    """Parse strict line-aware construction truth for proposed junctions."""
    with open_text(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(JUNCTION_TRUTH_COLUMNS):
            raise InputValidationError(f"{path}:1: unexpected junction truth columns")
        rows: list[ExpectedJunction] = []
        seen_cases: set[str] = set()
        seen_pairs: set[tuple[str, str]] = set()
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise InputValidationError(
                    f"{path}:{line_number}: row does not match junction truth columns"
                )
            try:
                parsed = _parse_junction_truth_row(row)
            except (KeyError, TypeError, ValueError) as error:
                raise InputValidationError(f"{path}:{line_number}: {error}") from error
            pair = (parsed.left_contig_id, parsed.right_contig_id)
            if parsed.case_id in seen_cases:
                raise InputValidationError(
                    f"{path}:{line_number}: duplicate junction case {parsed.case_id!r}"
                )
            if pair in seen_pairs:
                raise InputValidationError(
                    f"{path}:{line_number}: duplicate ordered junction pair {pair!r}"
                )
            seen_cases.add(parsed.case_id)
            seen_pairs.add(pair)
            rows.append(parsed)
    return tuple(sorted(rows, key=lambda item: item.case_id))


def evaluate_junction_support(
    evidence: tuple[TargetedJunctionEvidence, ...], policy: JunctionSupportPolicy
) -> JunctionSupportDecision:
    """Apply an explicit reviewed policy without authorizing any graph operation."""
    _validate_evidence_group(evidence)
    spanning = sum(item.spanning_reads for item in evidence)
    remapped = sum(len(item.remapped_read_names) for item in evidence)
    fraction = spanning / remapped if remapped else 0.0
    if not evidence:
        return JunctionSupportDecision(
            JunctionSupportStatus.DEFERRED,
            policy.technology,
            spanning,
            remapped,
            fraction,
            ("no targeted-remapping evidence supplied",),
        )
    if any(
        item.technology != policy.technology or item.remapping_preset != policy.remapping_preset
        for item in evidence
    ):
        return JunctionSupportDecision(
            JunctionSupportStatus.DEFERRED,
            policy.technology,
            spanning,
            remapped,
            fraction,
            ("evidence technology or remapping preset does not match the policy",),
        )
    if len(evidence) > 1:
        return JunctionSupportDecision(
            JunctionSupportStatus.DEFERRED,
            policy.technology,
            spanning,
            remapped,
            fraction,
            ("cross-sample junction evidence aggregation has not been reviewed",),
        )
    if not policy.reviewed:
        return JunctionSupportDecision(
            JunctionSupportStatus.DEFERRED,
            policy.technology,
            spanning,
            remapped,
            fraction,
            ("policy has not been reviewed against technology-specific junction truth",),
        )
    if evidence and any(
        item.minimum_spanning_flank < policy.minimum_spanning_flank for item in evidence
    ):
        return JunctionSupportDecision(
            JunctionSupportStatus.DEFERRED,
            policy.technology,
            spanning,
            remapped,
            fraction,
            ("one or more observations used a shorter spanning flank than policy requires",),
        )
    if evidence and any(
        item.minimum_mapping_quality != policy.minimum_mapping_quality for item in evidence
    ):
        return JunctionSupportDecision(
            JunctionSupportStatus.DEFERRED,
            policy.technology,
            spanning,
            remapped,
            fraction,
            ("evidence mapping-quality threshold does not match the policy",),
        )
    supported = (
        spanning >= policy.minimum_spanning_reads and fraction >= policy.minimum_spanning_fraction
    )
    return JunctionSupportDecision(
        JunctionSupportStatus.SUPPORTED if supported else JunctionSupportStatus.UNSUPPORTED,
        policy.technology,
        spanning,
        remapped,
        fraction,
        ("reviewed technology-specific thresholds satisfied",)
        if supported
        else ("reviewed technology-specific thresholds not satisfied",),
    )


def score_junction_observations(
    truth: tuple[ExpectedJunction, ...],
    observations: dict[str, tuple[TargetedJunctionEvidence, ...]],
) -> JunctionBenchmarkReport:
    """Score presence of spanning reads; this does not assess merge correctness."""
    known = {item.case_id for item in truth}
    unexpected = sorted(set(observations) - known)
    if unexpected:
        raise InputValidationError(f"junction observations contain unknown case: {unexpected[0]}")
    cases: list[JunctionBenchmarkCase] = []
    for expected in truth:
        reports = observations.get(expected.case_id, ())
        _validate_evidence_group(
            reports,
            expected_pair=(expected.left_contig_id, expected.right_contig_id),
        )
        spanning = sum(item.spanning_reads for item in reports)
        remapped = sum(len(item.remapped_read_names) for item in reports)
        detected = spanning > 0
        false_support = not expected.junction_is_true and detected
        missed = expected.junction_is_true and not detected
        absent = not reports
        reasons = [expected.reason]
        if absent:
            reasons.append("no targeted-remapping observation supplied")
        cases.append(
            JunctionBenchmarkCase(
                expected.case_id,
                expected.junction_is_true,
                spanning,
                remapped,
                detected,
                not false_support and not missed,
                false_support,
                missed,
                absent,
                tuple(reasons),
            )
        )
    ordered = tuple(sorted(cases, key=lambda item: item.case_id))
    summary = JunctionBenchmarkSummary(
        expected_cases=len(ordered),
        true_junctions=sum(item.junction_is_true for item in ordered),
        artificial_junctions=sum(not item.junction_is_true for item in ordered),
        observed_cases=sum(not item.absent_observation for item in ordered),
        absent_observations=sum(item.absent_observation for item in ordered),
        correct_detections=sum(item.correct_detection for item in ordered),
        false_support_cases=sum(item.false_support for item in ordered),
        missed_support_cases=sum(item.missed_support for item in ordered),
    )
    return JunctionBenchmarkReport(summary, ordered)


def _parse_junction_truth_row(row: dict[str, str]) -> ExpectedJunction:
    identifiers = (row["case_id"], row["left_contig"], row["right_contig"])
    if not all(identifiers):
        raise ValueError("junction case and contig identifiers cannot be empty")
    truth = _parse_boolean(row["junction_is_true"])
    circular = _parse_boolean(row["circular_wrap"])
    coordinate = int(row["source_coordinate"])
    spanning = int(row["selected_spanning_reads"])
    nonspanning = int(row["selected_nonspanning_reads"])
    if min(coordinate, spanning, nonspanning) < 0:
        raise ValueError("junction coordinates and read counts cannot be negative")
    if not row["expected_read_evidence"] or not row["reason"]:
        raise ValueError("junction evidence expectation and reason cannot be empty")
    return ExpectedJunction(
        identifiers[0],
        identifiers[1],
        identifiers[2],
        truth,
        coordinate,
        circular,
        spanning,
        nonspanning,
        row["expected_read_evidence"],
        row["reason"],
    )


def _parse_boolean(value: str) -> bool:
    if value not in {"true", "false"}:
        raise ValueError(f"invalid Boolean {value!r}")
    return value == "true"


def _validate_evidence_group(
    evidence: tuple[TargetedJunctionEvidence, ...],
    *,
    expected_pair: tuple[str, str] | None = None,
) -> None:
    if not evidence:
        return
    samples = [item.sample for item in evidence]
    if len(samples) != len(set(samples)):
        raise InputValidationError("junction evidence contains a duplicate sample")
    identities = {
        (
            item.left_contig_id,
            item.right_contig_id,
            item.provisional_reference_id,
            item.provisional_reference_length,
            item.provisional_reference_sha256,
            item.junction_position,
        )
        for item in evidence
    }
    if len(identities) != 1:
        raise InputValidationError(
            "junction evidence group does not describe one identical provisional junction"
        )
    pair = (evidence[0].left_contig_id, evidence[0].right_contig_id)
    if expected_pair is not None and pair != expected_pair:
        raise InputValidationError(
            f"junction observation pair {pair!r} does not match truth pair {expected_pair!r}"
        )
