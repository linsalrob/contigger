"""Transparent, domain-safe opening of plain and gzip-compressed text files."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

from contigger.exceptions import InputValidationError


@contextmanager
def open_text(
    path: Path, *, encoding: str = "utf-8", newline: str | None = None
) -> Iterator[TextIO]:
    """Open ``path`` as text, using gzip only for the explicit ``.gz`` suffix.

    Opening and deferred decompression errors are converted to stable domain errors.
    """
    try:
        if path.suffix.lower() == ".gz":
            with gzip.open(path, mode="rt", encoding=encoding, newline=newline) as handle:
                yield handle
        else:
            with path.open(mode="r", encoding=encoding, newline=newline) as handle:
                yield handle
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"cannot read text file {path}: {error}") from error
