"""Sweep candidate settings against checked-in Pseudomonas recall baselines."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from contigger.config import build_run_config
from contigger.pipeline_benchmark import evaluate_pipeline_benchmark


def main() -> int:
    """Evaluate deterministic candidate-setting combinations without changing baselines."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=Path("test_data"))
    parser.add_argument("--kmer-sizes", type=int, nargs="+", default=[21, 31])
    parser.add_argument("--window-sizes", type=int, nargs="+", default=[10, 15])
    parser.add_argument("--max-frequencies", type=int, nargs="+", default=[20, 50, 100])
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    config = build_run_config()
    rows: list[dict[str, object]] = []
    for kmer_size, window_size, maximum_frequency in itertools.product(
        arguments.kmer_sizes, arguments.window_sizes, arguments.max_frequencies
    ):
        reports = [
            evaluate_pipeline_benchmark(
                arguments.dataset,
                arguments.dataset / "alignments" / f"all_vs_all.{preset}.paf.gz",
                config,
                kmer_size=kmer_size,
                window_size=window_size,
                max_minimiser_frequency=maximum_frequency,
            )
            for preset in ("asm5", "asm20")
        ]
        rows.append(
            {
                "kmer_size": kmer_size,
                "window_size": window_size,
                "max_minimiser_frequency": maximum_frequency,
                "candidate_pairs": reports[0].summary.candidate_pairs,
                "candidate_recall_complete": all(
                    report.summary.candidate_stage_missed_relationships == 0 for report in reports
                ),
                "false_merges": [
                    report.summary.relationship_stage_false_merges for report in reports
                ],
                "missed_relationships": [
                    report.summary.relationship_stage_missed_relationships for report in reports
                ],
            }
        )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
