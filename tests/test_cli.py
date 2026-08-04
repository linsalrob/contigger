"""CLI contract tests."""

from pathlib import Path

import pytest

from contigger.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize("arguments", [["--help"], ["validate", "--help"], ["merge", "--help"]])
def test_help(arguments: list[str], capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as error:
        main(arguments)
    assert error.value.code == 0
    assert "usage:" in capsys.readouterr().out


def test_successful_manifest_validation(capsys: pytest.CaptureFixture[str]) -> None:
    status = main(["validate", "--manifest", str(FIXTURES / "samples.tsv")])
    assert status == 0
    assert "validated 2 sample(s)" in capsys.readouterr().out


def test_missing_manifest(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    status = main(["validate", "--manifest", str(tmp_path / "missing.tsv")])
    assert status != 0
    assert "cannot read manifest" in capsys.readouterr().err


def test_malformed_input(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    fasta = tmp_path / "bad.fa"
    fasta.write_text("sequence-before-header\n", encoding="utf-8")
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\tcontigs\nS1\tbad.fa\n", encoding="utf-8")
    assert main(["validate", "--manifest", str(manifest)]) != 0
    assert "before first header" in capsys.readouterr().err


def test_successful_merge_dry_run(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    prefix = tmp_path / "results" / "contigger"
    status = main(
        [
            "merge",
            "--manifest",
            str(FIXTURES / "samples.tsv"),
            "--output-prefix",
            str(prefix),
            "--dry-run",
        ]
    )
    captured = capsys.readouterr()
    assert status == 0
    assert '"mode": "dry-run"' in captured.out
    assert '"identity": 0.98' in captured.out
    assert captured.out.index('"sample": "S01"') < captured.out.index('"sample": "S02"')
    assert not prefix.parent.exists()


def test_real_merge_fails_clearly(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    status = main(
        [
            "merge",
            "--manifest",
            str(FIXTURES / "samples.tsv"),
            "--output-prefix",
            str(tmp_path / "contigger"),
        ]
    )
    assert status != 0
    assert "sequence merging is not implemented" in capsys.readouterr().err
