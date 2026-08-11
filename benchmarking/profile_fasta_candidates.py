"""Profile Contigger's catalogue and candidate stages on a real FASTA fixture."""

from __future__ import annotations

import argparse
import json
import platform
import resource
import time
from pathlib import Path

from contigger.catalogue import build_catalogue
from contigger.fasta import read_fasta
from contigger.minimisers import generate_candidates_with_metrics


def run(
    fasta: Path,
    *,
    sample: str,
    kmer_size: int,
    window_size: int,
    min_shared_minimisers: int,
    max_minimiser_frequency: int,
    terminal_band: int,
    candidate_shards: int,
    max_seed_pair_observations: int | None,
) -> dict[str, object]:
    """Load a real FASTA and measure the production catalogue/candidate implementation."""
    started = time.perf_counter()
    records = tuple(read_fasta(fasta, sample))
    load_seconds = time.perf_counter() - started
    started = time.perf_counter()
    catalogue = build_catalogue(records)
    catalogue_seconds = time.perf_counter() - started
    started = time.perf_counter()
    candidates, metrics = generate_candidates_with_metrics(
        catalogue.sequences,
        kmer_size=kmer_size,
        window_size=window_size,
        min_shared_minimisers=min_shared_minimisers,
        max_minimiser_frequency=max_minimiser_frequency,
        terminal_band=terminal_band,
        candidate_shards=candidate_shards,
        max_seed_pair_observations=max_seed_pair_observations,
    )
    candidate_seconds = time.perf_counter() - started
    lengths = sorted(record.length for record in records)
    input_count = len(records)
    return {
        "schema_version": 1,
        "source_fasta": str(fasta),
        "input_contigs": input_count,
        "input_bases": sum(record.length for record in records),
        "mean_contig_length": sum(lengths) / input_count,
        "median_contig_length": lengths[input_count // 2],
        "canonical_sequences": len(catalogue.sequences),
        "candidate_count": len(candidates),
        "candidate_reduction_factor": (
            input_count * (input_count - 1) / 2 / len(candidates) if candidates else None
        ),
        "timings_seconds": {
            "fasta_load": load_seconds,
            "catalogue": catalogue_seconds,
            "candidate_generation": candidate_seconds,
        },
        "candidate_configuration": {
            "kmer_size": kmer_size,
            "window_size": window_size,
            "min_shared_minimisers": min_shared_minimisers,
            "max_minimiser_frequency": max_minimiser_frequency,
            "terminal_band": terminal_band,
            "candidate_shards": candidate_shards,
            "max_seed_pair_observations": max_seed_pair_observations,
        },
        "candidate_generation": metrics.as_dict(),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "platform": platform.platform(),
    }


def main() -> int:
    """Profile a real fixture from command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--sample", default="real_scale")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kmer-size", type=int, default=21)
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--min-shared-minimisers", type=int, default=5)
    parser.add_argument("--max-minimiser-frequency", type=int, default=20)
    parser.add_argument("--terminal-band", type=int, default=1000)
    parser.add_argument("--candidate-shards", type=int, default=64)
    parser.add_argument("--max-seed-pair-observations", type=int)
    arguments = parser.parse_args()
    report = run(
        arguments.fasta,
        sample=arguments.sample,
        kmer_size=arguments.kmer_size,
        window_size=arguments.window_size,
        min_shared_minimisers=arguments.min_shared_minimisers,
        max_minimiser_frequency=arguments.max_minimiser_frequency,
        terminal_band=arguments.terminal_band,
        candidate_shards=arguments.candidate_shards,
        max_seed_pair_observations=arguments.max_seed_pair_observations,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
