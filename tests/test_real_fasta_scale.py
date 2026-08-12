"""Tests for deterministic, external real-FASTA scale fixture construction."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

from contigger.fasta import read_fasta

REPOSITORY_ROOT = Path(__file__).parents[1]
FIXTURE_SCRIPT = REPOSITORY_ROOT / "benchmarking" / "real_fasta_scale.py"
PROFILE_SCRIPT = REPOSITORY_ROOT / "benchmarking" / "profile_fasta_candidates.py"
MANIFEST_PROFILE_SCRIPT = REPOSITORY_ROOT / "benchmarking" / "profile_manifest_candidates.py"


def _run_fixture_script(source: Path, output: Path, *counts: int) -> dict[str, object]:
    """Run the fixture CLI as an external user would and return its metadata."""
    result = subprocess.run(
        [
            sys.executable,
            str(FIXTURE_SCRIPT),
            "--source",
            str(source),
            "--output-directory",
            str(output),
            "--records",
            *(str(count) for count in counts),
            "--seed",
            "test-seed",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(result.stdout)


def test_create_fixtures_are_nested_and_preserve_records(tmp_path: Path) -> None:
    """Hash-ranked fixtures are deterministic nested subsets with valid manifests."""
    source = tmp_path / "source.fasta"
    source.write_text(">a description\nAAAA\n>b\nCCCCCC\n>c\nGGG\n>d\nTTTTT\n", encoding="ascii")
    output = tmp_path / "fixtures"

    metadata = _run_fixture_script(source, output, 2, 3)

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
    """Every duplicate is rejected, including one outside the selected subset."""
    source = tmp_path / "duplicate.fasta"
    duplicate = "duplicate"
    unique = next(
        candidate
        for index in range(1000)
        if _priority(candidate := f"unique-{index}") < _priority(duplicate)
    )
    source.write_text(
        f">{duplicate}\nAAAA\n>{duplicate}\nCCCC\n>{unique}\nGGGG\n", encoding="ascii"
    )

    result = subprocess.run(
        [
            sys.executable,
            str(FIXTURE_SCRIPT),
            "--source",
            str(source),
            "--output-directory",
            str(tmp_path / "fixtures"),
            "--records",
            "1",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "duplicate FASTA identifier" in result.stderr


def test_create_fixtures_rejects_invalid_dna_symbols(tmp_path: Path) -> None:
    """A malformed source cannot produce a fixture that later fails Contigger validation."""
    source = tmp_path / "invalid.fasta"
    source.write_text(">invalid\nAAAACZCC\n", encoding="ascii")

    result = subprocess.run(
        [
            sys.executable,
            str(FIXTURE_SCRIPT),
            "--source",
            str(source),
            "--output-directory",
            str(tmp_path / "fixtures"),
            "--records",
            "1",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "invalid DNA symbol" in result.stderr


def test_profile_uses_real_fasta_records(tmp_path: Path) -> None:
    """The real-data profiler reports production candidate-stage counts."""
    fasta = tmp_path / "fixture.fasta"
    fasta.write_text(">a\nAAAACCCC\n>b\nGGGGTTTT\n", encoding="ascii")
    output = tmp_path / "profile.json"

    subprocess.run(
        [
            sys.executable,
            str(PROFILE_SCRIPT),
            "--fasta",
            str(fasta),
            "--output",
            str(output),
            "--kmer-size",
            "4",
            "--window-size",
            "2",
            "--min-shared-minimisers",
            "1",
            "--max-minimiser-frequency",
            "10",
            "--terminal-band",
            "10",
            "--candidate-shards",
            "2",
        ],
        check=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))

    assert result["input_contigs"] == 2
    assert result["input_bases"] == 16
    assert result["canonical_sequences"] == 1
    assert result["candidate_count"] == 0


def test_manifest_profile_validates_and_bounds_candidate_preflight(tmp_path: Path) -> None:
    """The manifest preflight uses production loading without invoking an aligner."""
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    first.write_text(">a\nAAAACCCC\n", encoding="ascii")
    second.write_text(">b\nGGGGTTTT\n", encoding="ascii")
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample\tcontigs\nfirst\tfirst.fasta\nsecond\tsecond.fasta\n",
        encoding="utf-8",
    )
    output = tmp_path / "manifest-profile.json"

    subprocess.run(
        [
            sys.executable,
            str(MANIFEST_PROFILE_SCRIPT),
            "--manifest",
            str(manifest),
            "--output",
            str(output),
            "--kmer-size",
            "4",
            "--window-size",
            "2",
            "--min-shared-minimisers",
            "1",
            "--max-minimiser-frequency",
            "10",
            "--terminal-band",
            "10",
            "--candidate-shards",
            "2",
            "--max-seed-pair-observations",
            "100",
            "--max-candidate-pairs",
            "1",
        ],
        check=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))

    assert result["input_samples"] == 2
    assert result["input_contigs"] == 2
    assert result["candidate_count"] == 0


def test_manifest_profile_handles_whitespace_only_fasta(tmp_path: Path) -> None:
    """A parser-accepted empty record set produces explicit nullable summaries."""
    fasta = tmp_path / "whitespace.fasta"
    fasta.write_text("\n \n", encoding="ascii")
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\tcontigs\nempty\twhitespace.fasta\n", encoding="utf-8")
    output = tmp_path / "empty-profile.json"

    subprocess.run(
        [
            sys.executable,
            str(MANIFEST_PROFILE_SCRIPT),
            "--manifest",
            str(manifest),
            "--output",
            str(output),
            "--max-seed-pair-observations",
            "1",
            "--max-candidate-pairs",
            "1",
        ],
        check=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["input_contigs"] == 0
    assert result["mean_contig_length"] is None
    assert result["median_contig_length"] is None


def _priority(identifier: str) -> int:
    """Return the fixture script's test-seed rank for a short test identifier."""
    payload = f"test-seed\0{identifier}".encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")
