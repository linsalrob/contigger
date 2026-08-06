"""CLI contract tests."""

from pathlib import Path

import pytest

from contigger.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize(
    "arguments",
    [
        ["--help"],
        ["validate", "--help"],
        ["merge", "--help"],
        ["classify-paf", "--help"],
        ["benchmark", "--help"],
        ["benchmark-pipeline", "--help"],
        ["catalogue", "--help"],
        ["candidates", "--help"],
    ],
)
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


def test_catalogue_writes_canonical_fasta_and_complete_provenance(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    fasta = tmp_path / "catalogue.fasta"
    provenance = tmp_path / "catalogue.provenance.tsv"
    status = main(
        [
            "catalogue",
            "--manifest",
            str(FIXTURES / "samples.tsv"),
            "--output-fasta",
            str(fasta),
            "--output-provenance",
            str(provenance),
        ]
    )
    assert status == 0
    assert fasta.read_text(encoding="ascii").startswith(">ctg_")
    assert len(provenance.read_text(encoding="utf-8").splitlines()) == 4
    assert "canonical sequence" in capsys.readouterr().out


def test_candidates_command_emits_evidence_not_relationships(tmp_path: Path) -> None:
    output = tmp_path / "candidates.tsv"
    status = main(
        [
            "candidates",
            "--manifest",
            str(FIXTURES / "samples.tsv"),
            "--output",
            str(output),
            "--kmer-size",
            "3",
            "--window-size",
            "2",
            "--min-shared-minimisers",
            "1",
            "--terminal-band",
            "4",
        ]
    )
    assert status == 0
    header = output.read_text(encoding="utf-8").splitlines()[0]
    assert "terminal_topologies" in header
    assert "relationship_type" not in header


def test_classify_paf_writes_deterministic_tsv(tmp_path: Path) -> None:
    paf = tmp_path / "hits.paf"
    paf.write_text(
        "z\t1000\t700\t1000\t+\ta\t1000\t0\t300\t300\t300\t60\n"
        "q\t1000\t300\t600\t+\tt\t1000\t400\t700\t300\t300\t255\n",
        encoding="utf-8",
    )
    output = tmp_path / "relationships.tsv"
    status = main(
        [
            "classify-paf",
            "--paf",
            str(paf),
            "--output",
            str(output),
            "--min-overlap",
            "100",
            "--min-containment",
            "50",
            "--end-tolerance",
            "10",
        ]
    )
    assert status == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[1].startswith("q\tt\tNO_RELATIONSHIP")
    assert lines[2].startswith("z\ta\tQUERY_SUFFIX_TO_TARGET_PREFIX")
    assert lines[2].split("\t")[13:15] == ["1", "0"]


def test_classify_paf_tsv_includes_accepted_and_rejected_reasons(tmp_path: Path) -> None:
    paf = tmp_path / "hits.paf"
    paf.write_text(
        "q\t1000\t700\t1000\t+\tt\t1000\t0\t300\t300\t300\t60\n"
        "q\t1000\t300\t600\t+\tt\t1000\t400\t700\t300\t300\t60\n",
        encoding="utf-8",
    )
    output = tmp_path / "relationships.tsv"
    status = main(
        [
            "classify-paf",
            "--paf",
            str(paf),
            "--output",
            str(output),
            "--min-overlap",
            "100",
            "--min-containment",
            "50",
            "--end-tolerance",
            "10",
        ]
    )
    assert status == 0
    reasons = output.read_text(encoding="utf-8").splitlines()[1].split("\t")[-1]
    assert reasons == "alignment lacks compatible terminal geometry; compatible terminal geometry"


def test_classify_paf_rejects_malformed_input(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    paf = tmp_path / "bad.paf"
    paf.write_text("bad\n", encoding="utf-8")
    status = main(["classify-paf", "--paf", str(paf), "--output", str(tmp_path / "out.tsv")])
    assert status != 0
    assert "PAF line 1" in capsys.readouterr().err


def test_classify_paf_rejects_impossible_alignment_spans(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    paf = tmp_path / "bad-spans.paf"
    paf.write_text("q\t0\t0\t0\t+\tt\t1000\t0\t300\t300\t300\t60\n", encoding="utf-8")
    status = main(["classify-paf", "--paf", str(paf), "--output", str(tmp_path / "out.tsv")])
    assert status != 0
    error = capsys.readouterr().err
    assert "PAF line 1" in error
    assert "matching bases" in error
