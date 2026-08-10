"""Checks for the portable Slurm accounting helper."""

import os
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts" / "collect_slurm_status.sh"


def test_slurm_status_script_has_valid_shell_syntax() -> None:
    subprocess.run(["bash", "-n", str(SCRIPT)], check=True)


def test_slurm_status_script_rejects_invalid_job_id(tmp_path: Path) -> None:
    result = subprocess.run(
        ["bash", str(SCRIPT), "not-a-job", str(tmp_path / "status.tsv")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "JOB_ID must be a numeric" in result.stderr


def test_slurm_status_uses_step_maxrss_and_records_sacct_provenance(tmp_path: Path) -> None:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_sacct = fake_bin / "sacct"
    fake_sacct.write_text(
        "#!/usr/bin/env bash\n"
        'if [[ "$1" == "--version" ]]; then echo \'slurm 24.05\'; exit 0; fi\n'
        'if [[ " $* " == *" -X "* ]]; then\n'
        "  echo '123|COMPLETED|0:0|00:10|32||58880M|cpu=32|cpu=32'\n"
        "else\n"
        "  echo '123|COMPLETED|0:0|00:10|32||58880M|cpu=32|cpu=32'\n"
        "  echo '123.batch|COMPLETED|0:0|00:10|32|64M|58880M|cpu=32|cpu=32'\n"
        "  echo '123.0|COMPLETED|0:0|00:10|32|128M|58880M|cpu=32|cpu=32'\n"
        "fi\n",
        encoding="utf-8",
    )
    fake_sacct.chmod(0o755)
    output = tmp_path / "status.tsv"
    environment = os.environ | {"PATH": f"{fake_bin}:{os.environ['PATH']}"}
    result = subprocess.run(
        ["bash", str(SCRIPT), "123", str(output)],
        capture_output=True,
        text=True,
        env=environment,
    )
    assert result.returncode == 0, result.stderr
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[5:12] == [
        "max_rss_kib",
        "max_rss",
        "req_mem",
        "req_tres",
        "alloc_tres",
        "sacct_version",
        "sacct_command",
    ]
    fields = lines[1].split("\t")
    assert fields[5:7] == ["131072", "128M"]
    assert fields[10] == "slurm 24.05"
    assert fields[11].startswith("sacct -X -j 123")
