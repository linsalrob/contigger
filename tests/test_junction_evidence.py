"""Targeted provisional-junction remapping tests without external tools."""

from __future__ import annotations

from pathlib import Path

import pytest

from contigger.evidence.bam import BamEvidenceProvider
from contigger.evidence.junctions import TargetedJunctionRemapper, _score_sam
from contigger.exceptions import InputValidationError
from contigger.models import (
    ContigEnd,
    JunctionRemappingRequest,
    SampleInput,
    SequenceRecord,
)
from contigger.utilities.subprocesses import CommandResult


class FakeTools:
    """Create deterministic samtools files and minimap2 SAM output."""

    def __init__(self) -> None:
        self.commands: list[tuple[str, ...]] = []

    def samtools(self, arguments: tuple[str, ...]) -> CommandResult:
        self.commands.append(arguments)
        operation = arguments[1]
        if operation == "--version":
            stdout = "samtools 1.22\n"
        elif arguments[1:3] == ("view", "-H"):
            stdout = "@SQ\tSN:left\tLN:100\n@SQ\tSN:right\tLN:100\n"
        elif operation == "view" and "-N" not in arguments:
            region = arguments[-1]
            name = "left-read" if region.startswith("left:") else "right-read"
            stdout = (
                f"{name}\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
                "shared\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
            )
        elif operation == "view" or operation == "collate":
            Path(arguments[arguments.index("-o") + 1]).touch()
            stdout = ""
        elif operation == "fastq":
            stdout = "@left-read\nACGT\n+\nIIII\n"
        else:
            stdout = ""
        return CommandResult(arguments, stdout, "", 0)

    def minimap2(self, arguments: tuple[str, ...]) -> CommandResult:
        self.commands.append(arguments)
        if arguments[1] == "--version":
            stdout = "2.30-r1287\n"
        else:
            stdout = (
                "@HD\tVN:1.6\n@SQ\tSN:provisional\tLN:180\n"
                "left-read\t0\tprovisional\t61\t60\t80M\t*\t0\t0\t*\t*\n"
                "right-read\t0\tprovisional\t91\t60\t60M\t*\t0\t0\t*\t*\n"
                "shared\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
            )
        return CommandResult(arguments, stdout, "", 0)


def setup_provider(tmp_path: Path, tools: FakeTools) -> BamEvidenceProvider:
    """Build a validated provider over placeholder alignment files."""
    fasta = tmp_path / "source.fasta"
    fasta.write_text(f">left\n{'A' * 100}\n>right\n{'C' * 100}\n", encoding="ascii")
    bam = tmp_path / "source.bam"
    bam.touch()
    Path(f"{bam}.bai").touch()
    executable = tmp_path / "samtools"
    executable.touch()
    return BamEvidenceProvider(
        SampleInput("sample-a", fasta, bam=bam, technology="ont"),
        executable=executable,
        runner=tools.samtools,
    )


def request() -> JunctionRemappingRequest:
    """Return one explicit provisional-reference request."""
    reference = SequenceRecord("provisional", "", "provisional", "", "A" * 180, 180)
    return JunctionRemappingRequest(
        "sample-a",
        "left",
        ContigEnd.SUFFIX,
        "right",
        ContigEnd.PREFIX,
        reference,
        100,
        extraction_distance=25,
        minimum_spanning_flank=20,
    )


def test_targeted_remapping_counts_distinct_flank_spanners(tmp_path: Path) -> None:
    tools = FakeTools()
    minimap = tmp_path / "minimap2"
    minimap.touch()
    evidence = TargetedJunctionRemapper(
        setup_provider(tmp_path, tools), minimap2=minimap, runner=tools.minimap2
    ).evaluate(request())
    assert evidence.selected_read_names == ("left-read", "right-read", "shared")
    assert evidence.technology == "ont"
    assert evidence.remapping_preset == "sr"
    assert evidence.remapped_read_names == ("left-read", "right-read")
    assert evidence.spanning_read_names == ("left-read",)
    assert evidence.spanning_reads == 1
    assert evidence.provisional_reference_length == 180
    assert evidence.provisional_reference_sha256 == (
        "dbcabec0bf27c25204ea59050e80eb4ffe6bc32a4e491338186fdcf62e49091e"
    )
    assert evidence.samtools_version == "samtools 1.22"
    assert evidence.minimap2_version == "2.30-r1287"
    assert "do not authorize a merge" in evidence.diagnostics[0]
    assert all(isinstance(command, tuple) for command in evidence.commands)


def test_repeated_evaluations_have_request_scoped_ordered_commands(tmp_path: Path) -> None:
    tools = FakeTools()
    minimap = tmp_path / "minimap2"
    minimap.touch()
    remapper = TargetedJunctionRemapper(
        setup_provider(tmp_path, tools), minimap2=minimap, runner=tools.minimap2
    )
    first = remapper.evaluate(request())
    second = remapper.evaluate(request())
    assert [Path(command[0]).name for command in second.commands] == [
        "samtools",
        "samtools",
        "minimap2",
        "samtools",
        "samtools",
        "samtools",
        "minimap2",
    ]
    assert len(first.commands) == len(second.commands) + 4
    assert [command[1] for command in second.commands] == [
        "view",
        "view",
        "--version",
        "view",
        "collate",
        "fastq",
        "-a",
    ]


def test_request_rejects_impossible_flank() -> None:
    reference = SequenceRecord("p", "", "p", "", "A" * 20, 20)
    with pytest.raises(InputValidationError, match="too short"):
        JunctionRemappingRequest(
            "s",
            "a",
            ContigEnd.SUFFIX,
            "b",
            ContigEnd.PREFIX,
            reference,
            10,
            minimum_spanning_flank=11,
        )


def test_sam_scoring_rejects_unselected_names_and_reference_mismatch() -> None:
    header = "@SQ\tSN:p\tLN:100\n"
    with pytest.raises(InputValidationError, match="unexpected read name"):
        _score_sam(
            header + "other\t0\tp\t1\t60\t50M\t*\t0\t0\t*\t*\n", "p", 100, 50, 10, ("chosen",)
        )
    with pytest.raises(InputValidationError, match="exact provisional reference"):
        _score_sam("@SQ\tSN:p\tLN:99\n", "p", 100, 50, 10, ())


def test_deletion_across_junction_is_not_spanning_support() -> None:
    sam = "@SQ\tSN:p\tLN:100\ngapped\t0\tp\t31\t60\t10M20D10M\t*\t0\t0\t*\t*\n"
    remapped, spanning = _score_sam(sam, "p", 100, 50, 5, ("gapped",))
    assert remapped == ("gapped",)
    assert spanning == ()
