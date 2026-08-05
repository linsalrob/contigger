"""Package-level regression tests for the checked-in Pseudomonas benchmark."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from contigger.benchmark import (
    BenchmarkReport,
    evaluate_benchmark,
    load_truth,
    write_json,
    write_tsv,
)
from contigger.cli import main
from contigger.config import build_run_config
from contigger.exceptions import InputValidationError

ROOT = Path(__file__).parents[1]
DATASET = ROOT / "test_data"
CONFIG = build_run_config()


def report(preset: str) -> BenchmarkReport:
    """Evaluate one checked-in PAF without invoking an external aligner."""
    paf = DATASET / "alignments" / f"all_vs_all.{preset}.paf.gz"
    return evaluate_benchmark(DATASET, paf, CONFIG)


def test_checked_in_manifest_with_gzip_fastas_validates() -> None:
    assert main(["validate", "--manifest", str(DATASET / "manifest.tsv")]) == 0


def test_truth_table_is_typed_sorted_and_complete() -> None:
    truth = load_truth(DATASET / "expected" / "expected_relationships.tsv")
    assert len(truth) == 74
    assert [(item.query_id, item.target_id) for item in truth] == sorted(
        (item.query_id, item.target_id) for item in truth
    )
    reverse = next(item for item in truth if item.case_id == "reverse_exact_r01")
    assert reverse.orientation is not None
    assert reverse.orientation.value == "-"


def test_truth_rejects_duplicate_ordered_pairs_with_line_number(tmp_path: Path) -> None:
    source = DATASET / "expected" / "expected_relationships.tsv"
    lines = source.read_text(encoding="utf-8").splitlines()
    path = tmp_path / "truth.tsv"
    path.write_text("\n".join([lines[0], lines[1], lines[1]]) + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match=r"truth\.tsv:3: duplicate ordered truth pair"):
        load_truth(path)


@pytest.mark.parametrize(
    ("preset", "correct", "false_merges", "missed", "boundary_failures", "self_hits"),
    [("asm5", 58, 2, 4, 6, 94), ("asm20", 60, 2, 2, 4, 92)],
)
def test_deterministic_checked_in_baseline(
    preset: str,
    correct: int,
    false_merges: int,
    missed: int,
    boundary_failures: int,
    self_hits: int,
) -> None:
    summary = report(preset).summary
    assert summary.correct_classifications == correct
    assert summary.false_merges == false_merges
    assert summary.missed_relationships == missed
    assert summary.threshold_boundary_failures == boundary_failures
    assert summary.self_alignments_excluded == self_hits
    assert summary.graph_level_ambiguity_cases_deferred == 4


def test_biological_regression_cases_are_explicit() -> None:
    cases = {case.case_id: case for case in report("asm20").cases if case.case_id}
    assert cases["circular_origin"].correct
    assert cases["internal_shared"].correct
    assert cases["reverse_exact_r03"].correct
    assert cases["repeat_ambiguity"].graph_level_deferred
    assert cases["identity_9800"].correct
    assert cases["length_1000"].correct
    assert cases["end_tolerance_51"].false_merge


def test_json_and_tsv_outputs_are_deterministic(tmp_path: Path) -> None:
    result = report("asm20")
    json_a, json_b = tmp_path / "a.json", tmp_path / "b.json"
    tsv_a, tsv_b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    outputs = ((json_a, write_json), (json_b, write_json), (tsv_a, write_tsv), (tsv_b, write_tsv))
    for path, writer in outputs:
        with path.open("w", encoding="utf-8", newline="") as output:
            writer(result, output)
    assert json_a.read_bytes() == json_b.read_bytes()
    assert tsv_a.read_bytes() == tsv_b.read_bytes()
    payload = json.loads(json_a.read_text(encoding="utf-8"))
    assert payload["summary"]["false_merges"] == 2


def test_fail_on_false_merge_exit_behaviour() -> None:
    arguments = [
        "benchmark",
        "--dataset",
        str(DATASET),
        "--paf",
        str(DATASET / "alignments" / "all_vs_all.asm20.paf.gz"),
        "--fail-on-false-merge",
    ]
    assert main(arguments) == 1


def test_fail_on_false_merge_succeeds_for_zero_count(monkeypatch: pytest.MonkeyPatch) -> None:
    real = report("asm20")
    zero = replace(real, summary=replace(real.summary, false_merges=0))
    monkeypatch.setattr("contigger.cli.evaluate_benchmark", lambda *_args: zero)
    assert main(["benchmark", "--dataset", ".", "--paf", "x", "--fail-on-false-merge"]) == 0
