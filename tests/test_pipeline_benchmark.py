"""Regression tests for the catalogue-to-relationship Pseudomonas baseline."""

from __future__ import annotations

import json
from io import StringIO
from pathlib import Path

import pytest

from contigger.cli import main
from contigger.config import build_run_config
from contigger.pipeline_benchmark import (
    PipelineBenchmarkReport,
    evaluate_pipeline_benchmark,
    write_pipeline_json,
)

DATASET = Path(__file__).parents[1] / "test_data"
BASELINE = DATASET.parent / "benchmarks" / "pseudomonas_pipeline_baseline.json"
CONFIG = build_run_config()


@pytest.fixture(scope="module")
def reports() -> dict[str, PipelineBenchmarkReport]:
    """Evaluate each checked-in PAF only once for this regression module."""
    return {
        preset: evaluate_pipeline_benchmark(
            DATASET, DATASET / "alignments" / f"all_vs_all.{preset}.paf.gz", CONFIG
        )
        for preset in ("asm5", "asm20")
    }


@pytest.mark.parametrize(
    ("preset", "missed"),
    (("asm5", 4), ("asm20", 2)),
)
def test_checked_in_pipeline_baseline(
    preset: str, missed: int, reports: dict[str, PipelineBenchmarkReport]
) -> None:
    report = reports[preset]
    summary = report.summary
    assert (summary.source_contigs, summary.canonical_sequences) == (90, 84)
    assert summary.candidate_pairs == 61
    assert summary.exact_truth_rows_recovered == summary.exact_truth_rows == 12
    assert summary.valid_nonexact_case_groups_recalled == 23
    assert summary.valid_nonexact_case_groups == 23
    assert summary.valid_nonexact_truth_rows_recalled == 40
    assert summary.candidate_stage_missed_relationships == 0
    assert summary.relationship_stage_false_merges == 2
    assert summary.relationship_stage_missed_relationships == missed
    assert summary.graph_level_ambiguity_cases_deferred == 4
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))[preset]
    assert baseline["paf_sha256"] == report.paf_checksum
    assert baseline["candidate_pairs"] == summary.candidate_pairs
    assert baseline["relationship_stage_false_merges"] == 2
    assert baseline["relationship_stage_missed_relationships"] == missed


def test_pipeline_json_is_deterministic(reports: dict[str, PipelineBenchmarkReport]) -> None:
    report = reports["asm5"]
    first = StringIO()
    second = StringIO()
    write_pipeline_json(report, first)
    write_pipeline_json(report, second)
    assert first.getvalue() == second.getvalue()
    assert "elapsed" not in first.getvalue()


def test_pipeline_cli_false_merge_policy_and_json(
    tmp_path: Path,
    reports: dict[str, PipelineBenchmarkReport],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "pipeline.json"
    monkeypatch.setattr(
        "contigger.cli.evaluate_pipeline_benchmark",
        lambda *_args, **_kwargs: reports["asm5"],
    )
    status = main(
        [
            "benchmark-pipeline",
            "--dataset",
            str(DATASET),
            "--paf",
            str(DATASET / "alignments" / "all_vs_all.asm5.paf.gz"),
            "--output-json",
            str(output),
            "--fail-on-false-merge",
        ]
    )
    assert status == 1
    assert (
        json.loads(output.read_text(encoding="utf-8"))["summary"]["relationship_stage_false_merges"]
        == 2
    )
