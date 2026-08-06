"""Sample-scoped BAM/CRAM validation and pileup tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from contigger.evidence.bam import BamEvidenceProvider
from contigger.exceptions import InputValidationError
from contigger.models import ContigEnd, SampleInput
from contigger.utilities.subprocesses import CommandResult


class FakeSamtools:
    """Deterministic command runner for samtools provider tests."""

    def __init__(self, header: str, pileup: str = "", view: str = "") -> None:
        self.header = header
        self.pileup_text = pileup
        self.view_text = view
        self.commands: list[tuple[str, ...]] = []

    def __call__(self, arguments: tuple[str, ...]) -> CommandResult:
        self.commands.append(arguments)
        if arguments[1] == "--version":
            stdout = "samtools 1.22\n"
        elif arguments[1:3] == ("view", "-H"):
            stdout = self.header
        elif arguments[1] == "mpileup":
            stdout = self.pileup_text
        elif arguments[1] == "view":
            if "-o" in arguments:
                Path(arguments[arguments.index("-o") + 1]).touch()
                stdout = ""
            else:
                stdout = self.view_text
        elif arguments[1] == "collate":
            Path(arguments[arguments.index("-o") + 1]).touch()
            stdout = ""
        elif arguments[1] == "fastq":
            stdout = "@read-a\nACGT\n+\nIIII\n"
        else:
            stdout = ""
        return CommandResult(arguments, stdout, "", 0)


def sample_files(tmp_path: Path, *, indexed: bool = True, suffix: str = ".bam") -> SampleInput:
    """Create source FASTA and placeholder BAM paths."""
    fasta = tmp_path / "source.fasta"
    fasta.write_text(">a\nACGT\n>b\nAAAAA\n", encoding="ascii")
    bam = tmp_path / f"source{suffix}"
    bam.touch()
    if indexed:
        Path(f"{bam}{'.bai' if suffix == '.bam' else '.crai'}").touch()
    return SampleInput("sample-a", fasta, bam=bam)


def provider(
    tmp_path: Path,
    runner: FakeSamtools,
    *,
    indexed: bool = True,
    suffix: str = ".bam",
) -> BamEvidenceProvider:
    """Build a provider with an injected executable and runner."""
    return BamEvidenceProvider(
        sample_files(tmp_path, indexed=indexed, suffix=suffix),
        executable=_fake_executable(tmp_path),
        runner=runner,
    )


def _fake_executable(tmp_path: Path) -> Path:
    executable = tmp_path / "samtools"
    executable.touch()
    return executable


def test_validate_source_checks_index_integrity_header_and_version(tmp_path: Path) -> None:
    runner = FakeSamtools("@HD\tVN:1.6\n@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n")
    evidence = provider(tmp_path, runner)
    assert evidence.validate_source() == ("a", "b")
    assert evidence.version == "samtools 1.22"
    assert evidence.commands == tuple(runner.commands)
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


def test_wrong_index_kind_does_not_satisfy_bam_requirement(tmp_path: Path) -> None:
    sample = sample_files(tmp_path, indexed=False)
    assert sample.bam is not None
    Path(f"{sample.bam}.crai").touch()
    evidence = BamEvidenceProvider(
        sample,
        executable=_fake_executable(tmp_path),
        runner=FakeSamtools("@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n"),
    )
    with pytest.raises(InputValidationError, match="lacks an adjacent index"):
        evidence.validate_source()


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


def test_zero_depth_placeholder_is_missing_evidence(tmp_path: Path) -> None:
    header = "@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n"
    rows = tuple(
        provider(tmp_path, FakeSamtools(header, "a\t1\tA\t0\t*\t*\t*\n")).pileup("a", 0, 1)
    )
    assert rows[0].depth == 0
    assert rows[0].allele_counts == {}
    assert rows[0].mean_base_quality is rows[0].mean_mapping_quality is None


def test_cram_commands_receive_sample_reference(tmp_path: Path) -> None:
    header = "@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n"
    evidence = provider(tmp_path, FakeSamtools(header), suffix=".cram")
    evidence.validate_source()
    tuple(evidence.pileup("a", 0, 1))
    view = next(command for command in evidence.commands if command[1] == "view")
    pileup = next(command for command in evidence.commands if command[1] == "mpileup")
    assert view[3:5] == ("-T", str(evidence.sample.contigs))
    assert "-f" in pileup
    assert str(evidence.sample.contigs) in pileup


def test_source_bam_never_claims_new_junction_support(tmp_path: Path) -> None:
    evidence = provider(tmp_path, FakeSamtools("@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n"))
    junction = evidence.junction_evidence("a", "b")
    assert not junction.testable
    assert "cannot validate" in junction.diagnostics[0]


def test_reads_near_end_are_primary_unique_and_deterministic(tmp_path: Path) -> None:
    sam = (
        "read-z\t0\ta\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "read-a\t99\ta\t1\t60\t4M\t=\t1\t0\tACGT\tIIII\n"
        "read-z\t0\ta\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
    )
    runner = FakeSamtools("@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n", view=sam)
    evidence = provider(tmp_path, runner)
    assert tuple(evidence.reads_near_end("b", ContigEnd.SUFFIX, 3)) == (
        "read-a",
        "read-z",
    )
    assert runner.commands[-1][-1] == "b:3-5"
    assert runner.commands[-1][runner.commands[-1].index("-F") + 1] == "2304"


def test_extract_reads_recovers_named_primary_records(tmp_path: Path) -> None:
    runner = FakeSamtools("@SQ\tSN:a\tLN:4\n@SQ\tSN:b\tLN:5\n")
    evidence = provider(tmp_path, runner)
    output = tmp_path / "selected.fastq"
    assert evidence.extract_reads(("read-z", "read-a", "read-z"), output) == (
        "read-a",
        "read-z",
    )
    assert output.is_file()
    selection = runner.commands[-3]
    assert ("-N" in selection) and ("-u" in selection) and ("-F" in selection)
    assert runner.commands[-2][1] == "collate"
