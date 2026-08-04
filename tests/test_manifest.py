"""Manifest parsing tests."""

from pathlib import Path

import pytest

from contigger.exceptions import ManifestError
from contigger.manifest import parse_manifest


def make_manifest(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "samples.tsv"
    path.write_text(body, encoding="utf-8")
    return path


def test_relative_paths_optional_fields_metadata_and_sorting(tmp_path: Path) -> None:
    (tmp_path / "a.fa").write_text(">a\nACGT\n", encoding="utf-8")
    (tmp_path / "reads.bam").write_bytes(b"placeholder")
    manifest = make_manifest(
        tmp_path,
        "sample\tcontigs\tbam\ttechnology\tassembly_graph\tgroup\n"
        " B \ta.fa\treads.bam\t illumina \t\t case \n"
        "A\ta.fa\t\t\t\tcontrol\n",
    )
    result = parse_manifest(manifest)
    assert [sample.sample for sample in result.samples] == ["A", "B"]
    assert result.samples[1].contigs == (tmp_path / "a.fa").resolve()
    assert result.samples[1].bam == (tmp_path / "reads.bam").resolve()
    assert result.samples[1].technology == "illumina"
    assert result.samples[1].metadata == {"group": "case"}
    assert "lacks an adjacent index" in result.warnings[0]


def test_missing_required_column_has_header_line(tmp_path: Path) -> None:
    with pytest.raises(ManifestError, match=r":1:.*contigs"):
        parse_manifest(make_manifest(tmp_path, "sample\nS1\n"))


def test_duplicate_sample_has_data_line(tmp_path: Path) -> None:
    manifest = make_manifest(tmp_path, "sample\tcontigs\nS1\ta.fa\nS1\tb.fa\n")
    with pytest.raises(ManifestError, match=r":3:.*duplicate"):
        parse_manifest(manifest, check_files=False)


def test_missing_required_value_has_data_line(tmp_path: Path) -> None:
    with pytest.raises(ManifestError, match=r":2:.*required"):
        parse_manifest(make_manifest(tmp_path, "sample\tcontigs\nS1\t\n"), check_files=False)
