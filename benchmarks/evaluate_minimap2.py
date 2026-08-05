#!/usr/bin/env python3
"""Evaluate minimap2 presets against deterministic conservative relationship truth."""

from __future__ import annotations

import argparse
import json
import tempfile
import time
from pathlib import Path
from typing import Any

from contigger.aligners.minimap2 import Minimap2Aligner
from contigger.config import build_run_config
from contigger.models import RelationshipType
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.synthetic import synthetic_cases
from contigger.utilities.subprocesses import find_executable


def evaluate(preset: str) -> dict[str, Any]:
    """Run one preset across independent fixture pairs and return benchmark counts."""
    executable = find_executable("minimap2")
    if executable is None:
        raise RuntimeError("minimap2 is not installed")
    counts: dict[str, Any] = {
        "preset": preset,
        "true_relationship_classifications": 0,
        "false_merges": 0,
        "missed_relationships": 0,
        "ambiguous_relationships": 0,
        "candidate_paf_records": 0,
        "elapsed_alignment_seconds": 0.0,
        "elapsed_classification_seconds": 0.0,
    }
    config = build_run_config(min_overlap=100, min_containment=100, end_tolerance=10)
    with tempfile.TemporaryDirectory(prefix="contigger-benchmark-") as directory:
        root = Path(directory)
        for case in synthetic_cases():
            target = root / "target.fa"
            query = root / "query.fa"
            target.write_text(f">target\n{case.target}\n", encoding="ascii")
            query.write_text(f">query\n{case.query}\n", encoding="ascii")
            aligner = Minimap2Aligner(executable=executable, preset=preset)
            started = time.perf_counter()
            hits = tuple(aligner.align_paths(target, query))
            counts["elapsed_alignment_seconds"] += time.perf_counter() - started
            counts["candidate_paf_records"] += len(hits)
            started = time.perf_counter()
            groups = list(group_ordered_pairs(hits))
            observed = (
                classify_pair(groups[0], config).relationship.relationship_type
                if groups
                else RelationshipType.NO_RELATIONSHIP
            )
            counts["elapsed_classification_seconds"] += time.perf_counter() - started
            if observed is case.truth:
                counts["true_relationship_classifications"] += 1
            if observed is RelationshipType.AMBIGUOUS_OVERLAP:
                counts["ambiguous_relationships"] += 1
            if case.truth in {
                RelationshipType.NO_RELATIONSHIP,
                RelationshipType.AMBIGUOUS_OVERLAP,
            } and observed not in {
                RelationshipType.NO_RELATIONSHIP,
                RelationshipType.AMBIGUOUS_OVERLAP,
            }:
                counts["false_merges"] += 1
            if (
                case.truth
                not in {
                    RelationshipType.NO_RELATIONSHIP,
                    RelationshipType.AMBIGUOUS_OVERLAP,
                }
                and observed is RelationshipType.NO_RELATIONSHIP
            ):
                counts["missed_relationships"] += 1
    return counts


def main() -> int:
    """Run selected presets and print false-merge-first human and optional JSON output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", action="append", choices=("asm5", "asm10", "asm20"))
    parser.add_argument("--json", type=Path, help="optional JSON report path")
    arguments = parser.parse_args()
    presets = arguments.preset or ["asm5", "asm20"]
    reports = [evaluate(preset) for preset in presets]
    for report in reports:
        print(f"preset: {report['preset']}")
        print(f"FALSE MERGES: {report['false_merges']}")
        for key in (
            "true_relationship_classifications",
            "missed_relationships",
            "ambiguous_relationships",
            "candidate_paf_records",
            "elapsed_alignment_seconds",
            "elapsed_classification_seconds",
        ):
            print(f"{key.replace('_', ' ')}: {report[key]}")
        print()
    if arguments.json:
        arguments.json.write_text(
            json.dumps(reports, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
