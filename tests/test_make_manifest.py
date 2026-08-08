"""Tests for the standalone manifest generator."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts" / "make_manifest.py"


def test_make_manifest_discovers_sidecars_and_relative_paths(tmp_path: Path) -> None:
    assemblies = tmp_path / "assemblies"
    assemblies.mkdir()
    (assemblies / "sample01.fasta.gz").write_text("compressed-placeholder", encoding="ascii")
    (assemblies / "sample01.gfa").write_text("H\tVN:Z:1.0\n", encoding="ascii")
    (assemblies / "sample01.bam").write_bytes(b"bam")
    (assemblies / "sample01.bam.bai").write_bytes(b"index")
    output = tmp_path / "manifests" / "samples.tsv"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "wrote 1 sample" in result.stdout
    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row == {
        "sample": "sample01",
        "contigs": str(assemblies / "sample01.fasta.gz"),
        "bam": str(assemblies / "sample01.bam"),
        "technology": "",
        "assembly_graph": str(assemblies / "sample01.gfa"),
    }


def test_make_manifest_rejects_colliding_sample_names(tmp_path: Path) -> None:
    (tmp_path / "sample.fa").write_text(">a\nACGT\n", encoding="ascii")
    (tmp_path / "sample.fasta").write_text(">a\nACGT\n", encoding="ascii")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "multiple FASTA files produce sample name" in result.stderr


def test_make_manifest_rejects_output_alias_and_normalizes_stems(tmp_path: Path) -> None:
    assembly = tmp_path / "sample.fa"
    assembly.write_text(">a\nACGT\n", encoding="ascii")
    alias = tmp_path / "manifest.tsv"
    alias.symlink_to(assembly)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--output", str(alias)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "aliases a discovered input" in result.stderr

    assembly.unlink()
    alias.unlink()
    (tmp_path / "sample .fa").write_text(">a\nACGT\n", encoding="ascii")
    output = tmp_path / "samples.tsv"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--output", str(output)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "sample\t" in output.read_text(encoding="utf-8")
