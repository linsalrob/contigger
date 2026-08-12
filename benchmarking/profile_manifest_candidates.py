"""Profile Contigger's catalogue and candidate stages for a complete manifest.

This deliberately stops before minimap2 alignment and graph construction.  It is for
bounded preflights of unfamiliar multi-assembly collections, not biological output.
"""

from __future__ import annotations

import argparse
import json
import platform
import resource
import time
from pathlib import Path

from contigger.catalogue import build_catalogue, load_source_sequences
from contigger.exceptions import InputValidationError
from contigger.manifest import parse_manifest
from contigger.minimisers import generate_candidates_with_metrics


def run(
    manifest: Path,
    *,
    kmer_size: int,
    window_size: int,
    min_shared_minimisers: int,
    max_minimiser_frequency: int,
    terminal_band: int,
    candidate_shards: int,
    max_seed_pair_observations: int,
    max_candidate_pairs: int,
) -> dict[str, object]:
    """Validate one manifest and measure the production candidate implementation."""
    if max_candidate_pairs < 1:
        raise InputValidationError("--max-candidate-pairs must be positive")
    started = time.perf_counter()
    validation = parse_manifest(manifest)
    records = load_source_sequences(validation.samples)
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
    if len(candidates) > max_candidate_pairs:
        raise InputValidationError(
            f"candidate pairs {len(candidates)} exceed --max-candidate-pairs "
            f"{max_candidate_pairs}; no alignment was started"
        )
    input_count = len(records)
    lengths = sorted(record.length for record in records)
    median_length: float | None
    if not lengths:
        median_length = None
    elif input_count % 2:
        median_length = float(lengths[input_count // 2])
    else:
        median_length = (lengths[input_count // 2 - 1] + lengths[input_count // 2]) / 2
    return {
        "schema_version": 1,
        "source_manifest": str(manifest.resolve()),
        "input_samples": len(validation.samples),
        "input_contigs": input_count,
        "input_bases": sum(record.length for record in records),
        "mean_contig_length": sum(lengths) / input_count if input_count else None,
        "median_contig_length": median_length,
        "canonical_sequences": len(catalogue.sequences),
        "candidate_count": len(candidates),
        "candidate_reduction_factor": (
            len(catalogue.sequences) * (len(catalogue.sequences) - 1) / 2 / len(candidates)
            if candidates
            else None
        ),
        "timings_seconds": {
            "manifest_and_fasta_load": load_seconds,
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
            "max_candidate_pairs": max_candidate_pairs,
        },
        "candidate_generation": metrics.as_dict(),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "platform": platform.platform(),
    }


def main() -> int:
    """Run a bounded candidate-only manifest preflight."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kmer-size", type=int, default=21)
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--min-shared-minimisers", type=int, default=5)
    parser.add_argument("--max-minimiser-frequency", type=int, default=20)
    parser.add_argument("--terminal-band", type=int, default=1000)
    parser.add_argument("--candidate-shards", type=int, default=64)
    parser.add_argument("--max-seed-pair-observations", type=int, required=True)
    parser.add_argument("--max-candidate-pairs", type=int, required=True)
    arguments = parser.parse_args()
    report = run(
        arguments.manifest,
        kmer_size=arguments.kmer_size,
        window_size=arguments.window_size,
        min_shared_minimisers=arguments.min_shared_minimisers,
        max_minimiser_frequency=arguments.max_minimiser_frequency,
        terminal_band=arguments.terminal_band,
        candidate_shards=arguments.candidate_shards,
        max_seed_pair_observations=arguments.max_seed_pair_observations,
        max_candidate_pairs=arguments.max_candidate_pairs,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
