"""Replaceable evidence-provider interface."""

from collections.abc import Iterable
from pathlib import Path
from typing import Protocol

from contigger.models import BaseEvidence, ContigEnd, JoinEvidence, SampleInput


class EvidenceProvider(Protocol):
    """Interface for sample-scoped source and junction evidence."""

    @property
    def sample(self) -> SampleInput:
        """Return the sample whose reads this provider represents."""
        ...

    def validate_source(self) -> tuple[str, ...]:
        """Validate alignment references against the source FASTA."""
        ...

    def pileup(self, contig_id: str, start: int, end: int) -> Iterable[BaseEvidence]:
        """Return evidence for a zero-based, half-open source interval."""
        ...

    def reads_near_end(self, contig_id: str, end: ContigEnd | str, distance: int) -> Iterable[str]:
        """Return stable identifiers for reads near a named contig end."""
        ...

    def extract_reads(self, read_names: Iterable[str], output_fastq: Path) -> tuple[str, ...]:
        """Recover selected reads and their primary mate records into FASTQ."""
        ...

    def junction_evidence(self, left_contig_id: str, right_contig_id: str) -> JoinEvidence:
        """Report whether targeted junction evidence is available."""
        ...
