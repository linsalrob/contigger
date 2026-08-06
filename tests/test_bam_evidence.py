"""Sample-scoped BAM/CRAM validation and pileup tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from contigger.evidence.bam import BamEvidenceProvider
from contigger.exceptions import FeatureNotImplementedError, InputValidationError
from contigger.models import SampleInput
from contigger.utilities.subprocesses import CommandResult


class FakeSamtools:
    """Deterministic command runner for samtools provider tests."""

    def __init__(self, header: str, pileup: str = "") -> None:
        self.header = header
        self.pileup_text = pileup
        self.commands: list[tuple[str, ...]] = []

    def __call__(self, arguments: tuple[str, ...]) -> CommandResult:
        self.commands.append(arguments)
        if arguments[1] == "--version":
            stdout = "samtools 1.22\n"
        elif arguments[1:3] == ("view", "-H"):
            stdout = self.header
        elif arguments[1] == "mpileup":
            stdout = self.pileup_text
        else:
            stdout = ""
        return CommandResult(arguments, stdout, "", 0)


def sample_files(tmp_path: Path, *, indexed: bool = True) -> SampleInput:
    """Create source FASTA and placeholder BAM paths."""
    fasta = tmp_path / "source.fasta"
    fasta.write_text(">a\nACGT\n>b\nAAAAA\n", encoding="ascii")
    bam = tmp_path / "source.bam"
    bam.touch()
    if indexed:
        Path(f"{bam}.bai").touch()
    return SampleInput("sample-a", fasta, bam=bam)


def provider(tmp_path: Path, runner: FakeSamtools, *, indexed: bool = True) -> BamEvidenceProvider:
    """Build a provider with an injected executable and runner."""
    return BamEvidenceProvider(
        sample_files(tmp_path, indexed=indexed),
        executable=Path("/test/samtools"),
        runner=runner,
    )


def test_validate_source_checks_index_integrity_header_and_version(tmp_path: Path) -> None:
    runner = FakeSamtools("@HD\tVN:1.6\n@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n")
    evidence = provider(tmp_path, runner)
    assert evidence.validate_source() == ("a", "b")
    assert evidence.version == "samtools 1.22"
    assert [command[1] for command in runner.commands] == [
        "--version",
        "quickcheck",
        "idxstats",
        "view",
    ]


def test_missing_index_and_reference_mismatch_are_errors(tmp_path: Path) -> None:
    runner = FakeSamtools("@SQ\tSN:a\tLN:4\n")
    with pytest.raises(InputValidationError, match="lacks an adjacent index"):
        provider(tmp_path, runner, indexed=False).validate_source()
    with pytest.raises(InputValidationError, match="missing reference 'b'"):
        provider(tmp_path, runner).validate_source()


def test_pileup_is_zero_based_sample_scoped_and_parses_markers(tmp_path: Path) -> None:
    runner = FakeSamtools(
        "@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n",
        "a\t2\tC\t3\t.,A\tABC\tDEF\na\t3\tG\t2\t^].+1t,$\tGH\tIJ\n",
    )
    evidence = provider(tmp_path, runner)
    rows = tuple(evidence.pileup("a", 1, 3))
    assert [(row.position, row.allele_counts, row.depth) for row in rows] == [
        (1, {"A": 1, "C": 2}, 3),
        (2, {"G": 2}, 2),
    ]
    assert all(row.sample == "sample-a" for row in rows)
    assert runner.commands[-1][1:] == (
        "mpileup",
        "-aa",
        "-s",
        "-r",
        "a:2-3",
        str(evidence.sample.bam),
    )


def test_invalid_interval_and_malformed_pileup_are_errors(tmp_path: Path) -> None:
    header = "@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n"
    with pytest.raises(InputValidationError, match="outside reference"):
        tuple(provider(tmp_path, FakeSamtools(header)).pileup("a", 0, 5))
    with pytest.raises(InputValidationError, match="truncated read start"):
        tuple(provider(tmp_path, FakeSamtools(header, "a\t1\tA\t1\t^\t!\t!\n")).pileup("a", 0, 1))


def test_source_bam_never_claims_new_junction_support(tmp_path: Path) -> None:
    evidence = provider(tmp_path, FakeSamtools("@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n"))
    junction = evidence.junction_evidence("a", "b")
    assert not junction.testable
    assert "cannot validate" in junction.diagnostics[0]
    with pytest.raises(FeatureNotImplementedError, match="read extraction"):
        tuple(evidence.reads_near_end("a", "suffix", 100))
