"""Checks for the portable Slurm accounting helper."""

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
