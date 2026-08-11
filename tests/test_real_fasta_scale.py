"""Tests for deterministic, external real-FASTA scale fixture construction."""

from __future__ import annotations

import json
from pathlib import Path

from benchmarking.profile_fasta_candidates import run
from benchmarking.real_fasta_scale import create_fixtures
from contigger.fasta import read_fasta


def test_create_fixtures_are_nested_and_preserve_records(tmp_path: Path) -> None:
    """Hash-ranked fixtures are deterministic nested subsets with valid manifests."""
    source = tmp_path / "source.fasta"
    source.write_text(">a description\nAAAA\n>b\nCCCCCC\n>c\nGGG\n>d\nTTTTT\n", encoding="ascii")
    output = tmp_path / "fixtures"

    metadata = create_fixtures(source, output, (2, 3), seed="test-seed")

    two = tuple(
        record.original_identifier for record in read_fasta(output / "real-contigs-2.fasta")
    )
    three = tuple(
        record.original_identifier for record in read_fasta(output / "real-contigs-3.fasta")
    )
    assert set(two) < set(three)
    assert metadata["source"] == {
        "path": str(source),
        "size_bytes": source.stat().st_size,
        "records": 4,
        "bases": 18,
    }
    assert (output / "real-contigs-2.samples.tsv").read_text(encoding="utf-8") == (
        "sample\tcontigs\nreal_contigs_2\treal-contigs-2.fasta\n"
    )
    assert (
        json.loads((output / "real-contigs-fixtures.json").read_text(encoding="utf-8"))["fixtures"][
            "3"
        ]["records"]
        == 3
    )


def test_create_fixtures_rejects_duplicate_identifiers(tmp_path: Path) -> None:
    """Duplicate source identifiers cannot yield a reproducible identifier-based fixture."""
    source = tmp_path / "duplicate.fasta"
    source.write_text(">same\nAAAA\n>same\nCCCC\n", encoding="ascii")

    try:
        create_fixtures(source, tmp_path / "fixtures", (1,), seed="test")
    except ValueError as error:
        assert "duplicate FASTA identifiers" in str(error)
    else:
        raise AssertionError("duplicate source identifiers were accepted")


def test_profile_uses_real_fasta_records(tmp_path: Path) -> None:
    """The real-data profiler reports production candidate-stage counts."""
    fasta = tmp_path / "fixture.fasta"
    fasta.write_text(">a\nAAAACCCC\n>b\nGGGGTTTT\n", encoding="ascii")

    result = run(
        fasta,
        sample="fixture",
        kmer_size=4,
        window_size=2,
        min_shared_minimisers=1,
        max_minimiser_frequency=10,
        terminal_band=10,
        candidate_shards=2,
        max_seed_pair_observations=None,
    )

    assert result["input_contigs"] == 2
    assert result["input_bases"] == 16
    assert result["canonical_sequences"] == 1
    assert result["candidate_count"] == 0
