"""Configuration normalisation at the public input boundary."""

from pathlib import Path

from contigger.exceptions import ConfigurationError
from contigger.models import ConflictPolicy, EvidenceMode, RunConfig


def normalise_percentage(value: float, name: str) -> float:
    """Normalise a human-readable percentage in [0, 100] to a fraction."""
    if not 0.0 <= value <= 100.0:
        raise ConfigurationError(f"{name} must be between 0 and 100")
    return value / 100.0


def build_run_config(
    *,
    identity: float = 98.0,
    min_overlap: int = 1000,
    min_containment: int = 500,
    containment_coverage: float = 98.0,
    end_tolerance: int = 50,
    kmer_size: int = 21,
    window_size: int = 10,
    min_shared_minimisers: int = 5,
    max_minimiser_frequency: int = 100,
    threads: int = 1,
    evidence: str = "none",
    conflict_policy: str = "reject",
    output_prefix: Path | str = Path("contigger"),
    deterministic_seed: int | None = None,
    emit_gfa: bool = False,
    minimap2_preset: str = "asm20",
    index_dir: Path | str | None = None,
    max_candidate_pairs: int | None = None,
    max_seed_pair_observations: int | None = None,
    candidate_shards: int = 16,
    max_queries_per_alignment_batch: int = 1000,
) -> RunConfig:
    """Validate raw CLI-style values and construct a normalised configuration."""
    try:
        evidence_mode = EvidenceMode(evidence)
    except ValueError as error:
        raise ConfigurationError(f"unsupported evidence mode: {evidence}") from error
    try:
        policy = ConflictPolicy(conflict_policy)
    except ValueError as error:
        raise ConfigurationError(f"unsupported conflict policy: {conflict_policy}") from error
    return RunConfig(
        identity=normalise_percentage(identity, "identity"),
        min_overlap=min_overlap,
        min_containment=min_containment,
        containment_coverage=normalise_percentage(containment_coverage, "containment coverage"),
        end_tolerance=end_tolerance,
        kmer_size=kmer_size,
        window_size=window_size,
        min_shared_minimisers=min_shared_minimisers,
        max_minimiser_frequency=max_minimiser_frequency,
        threads=threads,
        evidence=evidence_mode,
        conflict_policy=policy,
        output_prefix=Path(output_prefix),
        deterministic_seed=deterministic_seed,
        emit_gfa=emit_gfa,
        minimap2_preset=minimap2_preset,
        index_dir=Path(index_dir) if index_dir is not None else None,
        max_candidate_pairs=max_candidate_pairs,
        max_seed_pair_observations=max_seed_pair_observations,
        candidate_shards=candidate_shards,
        max_queries_per_alignment_batch=max_queries_per_alignment_batch,
    )
