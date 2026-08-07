"""Regression tests for the checked-in ONT junction-remapping baseline."""

from __future__ import annotations

import json
from pathlib import Path

from contigger.junction_remapping_benchmark import _build_requests

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "test_data"
BASELINE = ROOT / "benchmarks" / "pseudomonas_junction_remapping_baseline.json"


def test_provisional_junctions_match_construction_truth() -> None:
    requests = _build_requests(DATASET, 20)

    assert len(requests) == 19
    assert requests["circular_origin"].junction_position == 20_000
    assert requests["small_deletion"].provisional_reference.length == 35_000
    assert requests["small_insertion"].provisional_reference.length == 35_000


def test_checked_in_ont_remapping_baseline_is_conservative() -> None:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    summary = baseline["benchmark"]["summary"]

    assert baseline["dataset_version"] == "1.0.0"
    assert baseline["preset"] == "map-ont"
    assert baseline["minimum_spanning_flank"] == 20
    assert summary == {
        "absent_observations": 0,
        "artificial_junctions": 4,
        "correct_detections": 17,
        "expected_cases": 19,
        "false_support_cases": 0,
        "missed_support_cases": 2,
        "observed_cases": 19,
        "true_junctions": 15,
    }
    missed = {case["case_id"] for case in baseline["benchmark"]["cases"] if case["missed_support"]}
    assert missed == {"end_tolerance_49", "end_tolerance_50"}
