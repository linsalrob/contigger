"""Configuration boundary tests."""

import pytest

from contigger.config import build_run_config
from contigger.exceptions import ConfigurationError
from contigger.models import ConflictPolicy, EvidenceMode


def test_defaults_and_percentage_normalisation() -> None:
    config = build_run_config()
    assert config.identity == pytest.approx(0.98)
    assert config.containment_coverage == pytest.approx(0.98)
    assert config.kmer_size == 21
    assert config.evidence is EvidenceMode.NONE
    assert config.conflict_policy is ConflictPolicy.REJECT
    assert config.minimap2_preset == "asm20"
    assert config.index_dir is None


def test_index_and_preset_configuration() -> None:
    config = build_run_config(minimap2_preset="asm5", index_dir="cache")
    assert config.minimap2_preset == "asm5"
    assert config.index_dir.name == "cache"


def test_candidate_guard_configuration() -> None:
    config = build_run_config(max_candidate_pairs=12)
    assert config.max_candidate_pairs == 12
    with pytest.raises(ConfigurationError, match="candidate pairs"):
        build_run_config(max_candidate_pairs=0)


@pytest.mark.parametrize("identity", [-1.0, 100.1])
def test_invalid_identity(identity: float) -> None:
    with pytest.raises(ConfigurationError, match="identity"):
        build_run_config(identity=identity)


def test_invalid_thread_count() -> None:
    with pytest.raises(ConfigurationError, match="thread"):
        build_run_config(threads=0)


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("kmer_size", 0),
        ("window_size", 0),
        ("min_shared_minimisers", 0),
        ("max_minimiser_frequency", -1),
    ],
)
def test_invalid_minimiser_parameters(name: str, value: int) -> None:
    with pytest.raises(ConfigurationError):
        build_run_config(**{name: value})  # type: ignore[arg-type]


def test_unsupported_evidence_mode() -> None:
    with pytest.raises(ConfigurationError, match="evidence"):
        build_run_config(evidence="pooled")


def test_unsupported_conflict_policy() -> None:
    with pytest.raises(ConfigurationError, match="conflict"):
        build_run_config(conflict_policy="invent")
