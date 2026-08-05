import sys

import pytest
from validate_benchmark import resolve_tool


def test_validator_resolves_environment_override(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SAMTOOLS", sys.executable)
    assert resolve_tool("samtools") == sys.executable


def test_validator_reports_missing_tool_clearly(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SAMTOOLS", raising=False)
    monkeypatch.setenv("PATH", "")
    with pytest.raises(SystemExit, match="required installed tool not found on PATH: samtools"):
        resolve_tool("samtools")
