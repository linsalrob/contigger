"""Replaceable alignment interface."""

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Protocol

from contigger.models import AlignmentHit, SequenceRecord


class Aligner(Protocol):
    """Protocol implemented by native or external alignment engines."""

    @property
    def tool_name(self) -> str:
        """Return the alignment engine name."""
        ...

    @property
    def tool_version(self) -> str:
        """Return the exact detected tool version."""
        ...

    @property
    def last_command(self) -> tuple[str, ...] | None:
        """Return the last exact command for provenance."""
        ...

    def build_index(self, targets: Sequence[SequenceRecord], index_path: Path) -> Path:
        """Build or reuse an alignment index."""
        ...

    def align(
        self, queries: Iterable[SequenceRecord], targets: Iterable[SequenceRecord]
    ) -> Iterable[AlignmentHit]:
        """Yield typed alignments between supplied sequences."""
        ...


class IndexedAligner(Aligner, Protocol):
    """Alignment engine supporting an explicitly validated persistent index."""

    def align_indexed(
        self,
        queries: Iterable[SequenceRecord],
        targets: Sequence[SequenceRecord],
        index_path: Path,
    ) -> Iterable[AlignmentHit]:
        """Align queries to the exact targets represented by ``index_path``."""
        ...
