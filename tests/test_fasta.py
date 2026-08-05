"""FASTA syntax and sequence utility tests."""

import gzip
from pathlib import Path

import pytest

from contigger.exceptions import FastaFormatError
from contigger.fasta import read_fasta
from contigger.utilities.sequences import reverse_complement


def write_fasta(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "input.fasta"
    path.write_text(content, encoding="utf-8")
    return path


def test_multiline_description_and_lowercase(tmp_path: Path) -> None:
    records = list(read_fasta(write_fasta(tmp_path, ">id description words\nac gt\nnN\n"), "S1"))
    assert len(records) == 1
    assert records[0].original_identifier == "id"
    assert records[0].description == "description words"
    assert records[0].sequence == "ACGTNN"
    assert records[0].identifier == "S1:id"


@pytest.mark.parametrize(
    ("content", "message"),
    [
        (">same\nACG\n>same\nTTT\n", "duplicate"),
        (">empty\n>next\nACG\n", "empty sequence"),
        ("ACG\n>id\nACG\n", "before first header"),
        (">id\nAC-G\n", "invalid DNA"),
    ],
)
def test_rejects_malformed_fasta(tmp_path: Path, content: str, message: str) -> None:
    with pytest.raises(FastaFormatError, match=message):
        list(read_fasta(write_fasta(tmp_path, content)))


def test_reverse_complement_with_ambiguity() -> None:
    assert reverse_complement("ACGTRYSWKMBDHVN") == "NBDHVKMWSRYACGT"


def test_invalid_sequence_reports_last_sequence_line_before_blanks(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, ">bad\nAC-X\n\n\n>good\nACGT\n")

    with pytest.raises(FastaFormatError, match=r"input\.fasta:2: invalid DNA"):
        list(read_fasta(path))


def test_reads_gzip_fasta(tmp_path: Path) -> None:
    path = tmp_path / "input.fasta.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as output:
        output.write(">compressed\nacgt\n")
    assert [record.sequence for record in read_fasta(path)] == ["ACGT"]


def test_malformed_gzip_fasta_has_deterministic_domain_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.fasta.gz"
    path.write_bytes(b"not gzip")
    with pytest.raises(FastaFormatError, match=r"cannot read FASTA .*bad\.fasta\.gz"):
        list(read_fasta(path))
