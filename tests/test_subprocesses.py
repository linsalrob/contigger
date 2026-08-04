"""External command diagnostic tests."""

import sys

import pytest

from contigger.exceptions import ExternalToolError
from contigger.utilities.subprocesses import run_command


def test_failed_command_preserves_argument_vector() -> None:
    arguments = (
        sys.executable,
        "-c",
        "import sys; print('tool error', file=sys.stderr); raise SystemExit(7)",
        "argument with spaces",
        "semi;colon",
    )

    with pytest.raises(ExternalToolError) as error:
        run_command(arguments)

    message = str(error.value)
    assert "command failed (7)" in message
    assert repr(arguments) in message
    assert "tool error" in message
