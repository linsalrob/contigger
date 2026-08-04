"""BAM/CRAM provider scaffold with dependency-free path validation."""

from collections.abc import Iterable

from contigger.exceptions import FeatureNotImplementedError, InputValidationError
from contigger.models import BaseEvidence, JoinEvidence, SampleInput


class BamEvidenceProvider:
    """Future samtools/pysam-backed provider scoped to exactly one sample."""

    def __init__(self, sample: SampleInput) -> None:
        if sample.bam is None:
            raise InputValidationError(f"sample {sample.sample!r} has no BAM/CRAM input")
        if not sample.bam.is_file():
            raise InputValidationError(f"BAM/CRAM does not exist: {sample.bam}")
        self._sample = sample

    @property
    def sample(self) -> SampleInput:
        """Return the provider's sample input."""
        return self._sample

    def validate_source(self) -> tuple[str, ...]:
        """Reserve reference-name/length validation for a samtools-backed milestone."""
        raise FeatureNotImplementedError("BAM/CRAM reference validation is not implemented")

    def pileup(self, contig_id: str, start: int, end: int) -> Iterable[BaseEvidence]:
        """Reserve source-contig pileup analysis without inventing evidence."""
        raise FeatureNotImplementedError("BAM/CRAM pileup evidence is not implemented")

    def reads_near_end(self, contig_id: str, end: str, distance: int) -> Iterable[str]:
        """Reserve end-read extraction without inventing evidence."""
        raise FeatureNotImplementedError("read extraction near contig ends is not implemented")

    def junction_evidence(self, left_contig_id: str, right_contig_id: str) -> JoinEvidence:
        """State explicitly that a source BAM cannot yet validate a new junction."""
        return JoinEvidence(
            sample=self.sample.sample,
            left_contig_id=left_contig_id,
            right_contig_id=right_contig_id,
            testable=False,
            diagnostics=("targeted remapping is not implemented",),
        )
