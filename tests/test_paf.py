"""Strict PAF parsing tests that require no external aligner."""

import gzip
from io import StringIO
from pathlib import Path

import pytest

from contigger.aligners.minimap2 import parse_paf, parse_paf_line
from contigger.exceptions import InputValidationError
from contigger.models import AlignmentType, Orientation
from contigger.textio import open_text

PAF = "q\t1000\t700\t1000\t+\tt\t1000\t0\t300\t294\t300\t255"


def test_parse_paf_line_with_selected_optional_tags() -> None:
    hit = parse_paf_line(PAF + "\tAS:i:588\tcm:i:42\ts1:i:310\ts2:i:90\ttp:A:P\tzz:Z:kept")
    assert hit.orientation is Orientation.FORWARD
    assert hit.mapping_quality == 255
    assert hit.alignment_score == 588
    assert hit.supporting_seeds == 42
    assert hit.chaining_score == 310
    assert hit.secondary_chaining_score == 90
    assert hit.alignment_type is AlignmentType.PRIMARY


def test_parse_paf_stream_ignores_only_blank_lines() -> None:
    hits = list(parse_paf(StringIO("\n" + PAF + "\n  \n" + PAF.replace("q", "r", 1))))
    assert [hit.query_id for hit in hits] == ["q", "r"]


@pytest.mark.parametrize(
    "record",
    [
        "q\t100",
        PAF.replace("\t700\t", "\tnot-an-int\t", 1),
        PAF + "\tbroken",
        PAF + "\tAS:Z:wrong",
        PAF + "\ttp:A:X",
        PAF.rsplit("\t", 1)[0] + "\t256",
        PAF.replace("\t700\t1000\t", "\t1001\t1000\t", 1),
    ],
)
def test_malformed_paf_reports_physical_line_number(record: str) -> None:
    with pytest.raises(InputValidationError, match="PAF line 3"):
        list(parse_paf(["\n", "\n", record]))


def test_valid_paf_indel_can_have_block_longer_than_one_sequence_span() -> None:
    hit = parse_paf_line("q\t301\t0\t301\t+\tt\t300\t0\t300\t300\t301\t60")
    assert hit.alignment_block_length == 301


@pytest.mark.parametrize("compressed", [False, True])
def test_parse_plain_and_gzip_paf(tmp_path: Path, compressed: bool) -> None:
    path = tmp_path / ("hits.paf.gz" if compressed else "hits.paf")
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as output:
            output.write(PAF + "\n")
    else:
        path.write_text(PAF + "\n", encoding="utf-8")
    with open_text(path) as handle:
        assert [hit.query_id for hit in parse_paf(handle)] == ["q"]


def test_malformed_gzip_paf_has_informative_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.paf.gz"
    path.write_bytes(b"not gzip")
    with (
        pytest.raises(InputValidationError, match=r"cannot read text file .*bad\.paf\.gz"),
        open_text(path) as handle,
    ):
        list(parse_paf(handle))
