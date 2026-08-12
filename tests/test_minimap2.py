"""Minimap2 adapter unit tests using mocked command execution."""

from pathlib import Path

import pytest

import contigger.aligners.minimap2 as minimap2_module
from contigger.aligners.minimap2 import Minimap2Aligner
from contigger.exceptions import InputValidationError
from contigger.models import AlignmentHit, Orientation, SequenceRecord
from contigger.utilities.subprocesses import CommandResult


def test_align_paths_captures_version_command_and_parses_paf(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    target = tmp_path / "target.fa"
    query = tmp_path / "query.fa"
    target.write_text(">t\nAAAA\n", encoding="ascii")
    query.write_text(">q\nAAAA\n", encoding="ascii")
    calls: list[tuple[str, ...]] = []

    def fake_run(command: tuple[str, ...]) -> CommandResult:
        calls.append(command)
        stdout = (
            "2.28-r1209\n"
            if command[-1] == "--version"
            else "q\t4\t0\t4\t+\tt\t4\t0\t4\t4\t4\t255\ttp:A:P\n"
        )
        return CommandResult(command, stdout, "", 0)

    monkeypatch.setattr(minimap2_module, "run_command", fake_run)
    aligner = Minimap2Aligner(Path("/opt/minimap2"), threads=3, preset="asm5")
    hits = list(aligner.align_paths(target, query))
    assert len(hits) == 1
    assert aligner.tool_version == "2.28-r1209"
    assert calls[1] == ("/opt/minimap2", "-x", "asm5", "-t", "3", str(target), str(query))
    assert aligner.last_command == calls[1]


@pytest.mark.parametrize("preset", ["map-ont", "unknown"])
def test_rejects_non_assembly_presets(preset: str) -> None:
    with pytest.raises(InputValidationError, match="preset"):
        Minimap2Aligner(Path("minimap2"), preset=preset)


def test_rejects_invalid_threads() -> None:
    with pytest.raises(InputValidationError, match="thread"):
        Minimap2Aligner(Path("minimap2"), threads=0)


def test_typed_alignment_materialises_safe_temporary_fastas(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    commands: list[tuple[str, ...]] = []

    def fake_run(command: tuple[str, ...]) -> CommandResult:
        commands.append(command)
        stdout = (
            "2.31-r1302\n"
            if command[-1] == "--version"
            else "q\t4\t0\t4\t+\tt\t4\t0\t4\t4\t4\t60\n"
        )
        return CommandResult(command, stdout, "", 0)

    monkeypatch.setattr(minimap2_module, "run_command", fake_run)
    query = SequenceRecord("q", "", "q", "", "ACGT", 4)
    target = SequenceRecord("t", "", "t", "", "ACGT", 4)
    aligner = Minimap2Aligner(Path("/opt/minimap2"), preset="asm5")
    hits = tuple(aligner.align((query,), (target,)))
    assert [(hit.query_id, hit.target_id) for hit in hits] == [("q", "t")]
    assert commands[1][1:3] == ("-x", "asm5")
    assert commands[1][-2].endswith("targets.fasta")
    assert commands[1][-1].endswith("queries.fasta")


def test_build_index_reuses_only_matching_content(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls: list[tuple[str, ...]] = []

    def fake_run(command: tuple[str, ...]) -> CommandResult:
        calls.append(command)
        if command[-1] == "--version":
            return CommandResult(command, "2.31-r1302\n", "", 0)
        Path(command[command.index("-d") + 1]).write_bytes(b"mmi")
        return CommandResult(command, "", "", 0)

    monkeypatch.setattr(minimap2_module, "run_command", fake_run)
    aligner = Minimap2Aligner(Path("/opt/minimap2"))
    index = tmp_path / "target.mmi"
    target = SequenceRecord("t", "", "t", "", "ACGT", 4)
    assert aligner.build_index((target,), index) == index
    assert aligner.build_index((target,), index) == index
    reused = Minimap2Aligner(Path("/opt/minimap2"))
    assert reused.build_index((target,), index) == index
    assert reused.build_index((target,), index) == index
    assert reused.index_reuses == 1
    assert len(calls) == 3
    assert calls[1][1:4] == ("-x", "asm20", "-d")
    with pytest.raises(InputValidationError, match="does not match"):
        aligner.build_index((SequenceRecord("t", "", "t", "", "AAAA", 4),), index)
    with pytest.raises(InputValidationError, match="does not match"):
        Minimap2Aligner(Path("/opt/minimap2"), preset="asm5").build_index((target,), index)
    index.write_bytes(b"")
    with pytest.raises(InputValidationError, match="empty"):
        aligner.build_index((target,), index)


def test_incomplete_index_is_rejected(tmp_path: Path) -> None:
    index = tmp_path / "target.mmi"
    index.write_bytes(b"mmi")
    aligner = Minimap2Aligner(Path("/opt/minimap2"))
    target = SequenceRecord("t", "", "t", "", "ACGT", 4)
    with pytest.raises(InputValidationError, match="incomplete"):
        aligner.build_index((target,), index)


def test_metadata_failure_does_not_publish_index(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def fake_run(command: tuple[str, ...]) -> CommandResult:
        if command[-1] == "--version":
            return CommandResult(command, "2.31-r1302\n", "", 0)
        Path(command[command.index("-d") + 1]).write_bytes(b"mmi")
        return CommandResult(command, "", "", 0)

    original_write_text = Path.write_text

    def fail_metadata(path: Path, data: str, **kwargs: object) -> int:
        if path.name == "target.mmi.json":
            raise OSError("simulated metadata failure")
        return original_write_text(path, data, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(minimap2_module, "run_command", fake_run)
    monkeypatch.setattr(Path, "write_text", fail_metadata)
    index = tmp_path / "target.mmi"
    target = SequenceRecord("t", "", "t", "", "ACGT", 4)
    with pytest.raises(InputValidationError, match="simulated metadata failure"):
        Minimap2Aligner(Path("/opt/minimap2")).build_index((target,), index)
    assert not index.exists()
    assert not index.with_name("target.mmi.json").exists()


def test_align_indexed_rejects_unexpected_query_identifier(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    query = SequenceRecord("q", "", "q", "", "ACGT", 4)
    target = SequenceRecord("t", "", "t", "", "ACGT", 4)
    hit = AlignmentHit(
        "unexpected",
        "t",
        4,
        4,
        0,
        4,
        0,
        4,
        Orientation.FORWARD,
        4,
        4,
    )
    monkeypatch.setattr(Minimap2Aligner, "build_index", lambda *_args: tmp_path / "x.mmi")
    monkeypatch.setattr(Minimap2Aligner, "align_paths", lambda *_args: iter((hit,)))
    aligner = Minimap2Aligner(Path("/opt/minimap2"))
    with pytest.raises(InputValidationError, match="unexpected query"):
        tuple(aligner.align_indexed((query,), (target,), tmp_path / "x.mmi"))
