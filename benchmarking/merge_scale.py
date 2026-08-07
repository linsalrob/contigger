"""Optional deterministic catalogue/candidate scale benchmark.

This is intentionally outside ordinary unit-test execution.  It reports the
pre-alignment stages so larger runs can be profiled without committing data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from contigger.catalogue import build_catalogue
from contigger.minimisers import generate_candidates
from contigger.models import SequenceRecord


def run(contig_count: int, *, sequence_length: int = 1000) -> dict[str, object]:
    """Generate deterministic contigs and measure catalogue/minimiser stages."""
    if contig_count < 1 or sequence_length < 21:
        raise ValueError("contig count must be positive and sequence length at least 21")
    records = tuple(
        SequenceRecord(
            f"S{i % 10:02d}:synthetic_{i:07d}",
            f"S{i % 10:02d}",
            f"synthetic_{i:07d}",
            "deterministic scale benchmark",
            "".join(
                "ACGT"[value % 4]
                for value in hashlib.sha256(f"contig-{i}".encode()).digest()
                for _ in range((sequence_length + 31) // 32)
            )[:sequence_length],
            sequence_length,
            i,
        )
        for i in range(contig_count)
    )
    started = time.perf_counter()
    catalogue = build_catalogue(records)
    catalogue_seconds = time.perf_counter() - started
    started = time.perf_counter()
    candidates = generate_candidates(
        catalogue.sequences,
        kmer_size=21,
        window_size=10,
        min_shared_minimisers=5,
        max_minimiser_frequency=100,
        terminal_band=1000,
    )
    minimiser_seconds = time.perf_counter() - started
    return {
        "input_contigs": contig_count,
        "total_bases": contig_count * sequence_length,
        "canonical_sequences": len(catalogue.sequences),
        "catalogue_seconds": catalogue_seconds,
        "minimiser_seconds": minimiser_seconds,
        "candidate_count": len(candidates),
        "candidate_reduction_factor": (
            (contig_count * (contig_count - 1) / 2) / len(candidates) if candidates else None
        ),
    }


def main() -> int:
    """Run the optional scale benchmark and write JSON."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--contigs", type=int, default=10_000)
    parser.add_argument("--sequence-length", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = run(arguments.contigs, sequence_length=arguments.sequence_length)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
