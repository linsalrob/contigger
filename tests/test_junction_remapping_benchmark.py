"""Regression tests for the checked-in ONT junction-remapping baseline."""

from __future__ import annotations

import json
from pathlib import Path

from contigger.junction_remapping_benchmark import (
    JunctionRemappingBenchmarkReport,
    SampleScopedJunctionCase,
    _build_requests,
    benchmark_junction_policy_candidates,
)

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "test_data"
BASELINE = ROOT / "benchmarks" / "pseudomonas_junction_remapping_baseline.json"
CONFIGURATION_BASELINE = ROOT / "benchmarks" / "pseudomonas_junction_configuration_baseline.json"
POLICY_BASELINE = ROOT / "benchmarks" / "pseudomonas_junction_policy_candidate_baseline.json"


def test_provisional_junctions_match_construction_truth() -> None:
    requests = _build_requests(DATASET, 20)

    assert len(requests) == 19
    assert requests["circular_origin"].junction_position == 20_000
    assert requests["small_deletion"].provisional_reference.length == 35_000
    assert requests["small_insertion"].provisional_reference.length == 35_000


def test_checked_in_ont_remapping_baseline_is_conservative() -> None:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    summary = baseline["summary"]

    assert baseline["dataset_version"] == "1.0.0"
    assert baseline["preset"] == "map-ont"
    assert baseline["minimum_spanning_flank"] == 20
    assert baseline["minimum_mapping_quality"] == 0
    assert summary == {
        "artificial_sample_cases": 4,
        "expected_junctions": 19,
        "expected_sample_cases": 49,
        "false_support_baseline_established": True,
        "false_support_sample_cases": 0,
        "missed_support_sample_cases": 6,
        "testable_artificial_sample_cases": 4,
        "true_sample_cases": 45,
    }
    assert baseline["aggregate_case_status"]["missed_true_cases"] == [
        "end_tolerance_49",
        "end_tolerance_50",
    ]
    assert baseline["aggregate_case_status"]["untestable_artificial_cases"] == []
    controls = [
        case
        for case in baseline["cases"]
        if case["observation_kind"] == "synthetic-source-end-control"
    ]
    assert len(controls) == 4
    assert all(case["remapped_reads"] > 0 and case["spanning_reads"] == 0 for case in controls)


def test_checked_in_configuration_matrix_remains_unreviewed() -> None:
    baseline = json.loads(CONFIGURATION_BASELINE.read_text(encoding="utf-8"))

    assert baseline["configurations"] == [
        {"minimum_mapping_quality": 0, "minimum_spanning_flank": 20},
        {"minimum_mapping_quality": 20, "minimum_spanning_flank": 20},
        {"minimum_mapping_quality": 0, "minimum_spanning_flank": 100},
        {"minimum_mapping_quality": 20, "minimum_spanning_flank": 100},
    ]
    assert baseline["false_support_cases_each"] == 0
    assert baseline["testable_artificial_controls_each"] == 4
    assert baseline["policy_reviewed"] is False


def test_policy_candidate_benchmark_preserves_negative_controls() -> None:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    cases = tuple(SampleScopedJunctionCase(**case) for case in baseline["cases"])
    report = JunctionRemappingBenchmarkReport(
        baseline["dataset_version"],
        baseline["dataset_truth_sha256"],
        baseline["preset"],
        baseline["minimum_spanning_flank"],
        baseline["minimum_mapping_quality"],
        baseline["minimap2_version"],
        baseline["contigger_version"],
        None,  # type: ignore[arg-type]
        cases,
        baseline["aggregate_case_status"],
        baseline["observations"],
    )
    results = benchmark_junction_policy_candidates(report, ((1, 0.0), (3, 0.3), (5, 0.5)))

    assert [(item.false_support_cases, item.missed_true_cases) for item in results] == [
        (0, 6),
        (0, 9),
        (0, 12),
    ]
    assert all(not item.reviewed for item in results)
    policy_baseline = json.loads(POLICY_BASELINE.read_text(encoding="utf-8"))
    assert [
        {
            "minimum_spanning_fraction": item.minimum_spanning_fraction,
            "minimum_spanning_reads": item.minimum_spanning_reads,
            "false_support_cases": item.false_support_cases,
            "missed_true_cases": item.missed_true_cases,
            "supported_true_cases": item.supported_true_cases,
        }
        for item in results
    ] == policy_baseline["candidates"]
    assert policy_baseline["negative_controls"] == 4
    assert policy_baseline["reviewed"] is False
