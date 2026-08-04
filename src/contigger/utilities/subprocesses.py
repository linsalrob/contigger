"""Safe external command execution and provenance capture."""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from contigger.exceptions import ExternalToolError


@dataclass(frozen=True, slots=True)
class CommandResult:
    """Captured result of an external command."""

    command: tuple[str, ...]
    stdout: str
    stderr: str
    returncode: int


def find_executable(name: str) -> Path | None:
    """Return an absolute executable path when available."""
    found = shutil.which(name)
    return Path(found) if found else None


def run_command(arguments: Sequence[str]) -> CommandResult:
    """Run an argument array without shell interpolation and capture all output."""
    if not arguments:
        raise ExternalToolError("external command cannot be empty")
    try:
        completed = subprocess.run(  # noqa: S603 - validated argument array is intentional
            list(arguments), check=False, capture_output=True, text=True
        )
    except OSError as error:
        raise ExternalToolError(f"could not execute {arguments[0]!r}: {error}") from error
    result = CommandResult(
        tuple(arguments), completed.stdout, completed.stderr, completed.returncode
    )
    if completed.returncode != 0:
        raise ExternalToolError(
            f"command failed ({completed.returncode}): {' '.join(arguments)}\n"
            f"{completed.stderr.strip()}"
        )
    return result
