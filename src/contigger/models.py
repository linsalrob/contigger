"""Stable typed data models shared across Contigger components."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any

from contigger.exceptions import ConfigurationError, InputValidationError


class Orientation(StrEnum):
    """Explicit orientation of a query relative to a target."""

    FORWARD = "+"
    REVERSE = "-"


class RelationshipType(StrEnum):
    """Supported pairwise sequence relationship classifications."""

    EXACT_MATCH = "EXACT_MATCH"
    QUERY_CONTAINED_IN_TARGET = "QUERY_CONTAINED_IN_TARGET"
    TARGET_CONTAINED_IN_QUERY = "TARGET_CONTAINED_IN_QUERY"
    QUERY_SUFFIX_TO_TARGET_PREFIX = "QUERY_SUFFIX_TO_TARGET_PREFIX"
    TARGET_SUFFIX_TO_QUERY_PREFIX = "TARGET_SUFFIX_TO_QUERY_PREFIX"
    AMBIGUOUS_OVERLAP = "AMBIGUOUS_OVERLAP"
    NO_RELATIONSHIP = "NO_RELATIONSHIP"


class EvidenceMode(StrEnum):
    """Evidence source made available to a future decision policy."""

    NONE = "none"
    ALIGNMENTS = "alignments"
    READS = "reads"


class ConflictPolicy(StrEnum):
    """Declared policy for future sequence conflicts."""

    REPRESENTATIVE = "representative"
    MAJORITY = "majority"
    QUALITY_WEIGHTED = "quality-weighted"
    SAMPLE_AWARE = "sample-aware"
    AMBIGUOUS = "ambiguous"
    REJECT = "reject"


@dataclass(frozen=True, slots=True)
class SampleInput:
    """Files and metadata belonging to one biological sample."""

    sample: str
    contigs: Path
    bam: Path | None = None
    technology: str | None = None
    assembly_graph: Path | None = None
    metadata: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """A validated source sequence with stable source identity."""

    identifier: str
    source_sample: str
    original_identifier: str
    description: str
    sequence: str
    length: int
    source_ordinal: int | None = None

    def __post_init__(self) -> None:
        if self.length != len(self.sequence):
            raise InputValidationError(
                f"declared length {self.length} does not match sequence length {len(self.sequence)}"
            )


@dataclass(frozen=True, slots=True)
class CandidatePair:
    """A positional seed-supported request for pairwise alignment."""

    query_id: str
    target_id: str
    shared_minimisers: int
    orientation: Orientation | None = None


@dataclass(frozen=True, slots=True)
class AlignmentHit:
    """A zero-based, half-open pairwise alignment observation."""

    query_id: str
    target_id: str
    query_length: int
    target_length: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    orientation: Orientation
    matching_bases: int
    alignment_block_length: int
    mapping_quality: int | None = None
    alignment_score: int | None = None
    supporting_seeds: int | None = None

    def __post_init__(self) -> None:
        dimensions = (self.query_length, self.target_length, self.alignment_block_length)
        if any(value < 0 for value in dimensions):
            raise InputValidationError("alignment lengths cannot be negative")
        if not 0 <= self.query_start <= self.query_end <= self.query_length:
            raise InputValidationError("query coordinates are outside the query sequence")
        if not 0 <= self.target_start <= self.target_end <= self.target_length:
            raise InputValidationError("target coordinates are outside the target sequence")
        if not 0 <= self.matching_bases <= self.alignment_block_length:
            raise InputValidationError("matching bases must be within the alignment block")

    @property
    def identity(self) -> float:
        """Return matching bases divided by alignment block length, or zero."""
        if not self.alignment_block_length:
            return 0.0
        return self.matching_bases / self.alignment_block_length

    @property
    def query_coverage(self) -> float:
        """Return the aligned query span as a fraction of query length."""
        return (self.query_end - self.query_start) / self.query_length if self.query_length else 0.0

    @property
    def target_coverage(self) -> float:
        """Return the aligned target span as a fraction of target length."""
        if not self.target_length:
            return 0.0
        return (self.target_end - self.target_start) / self.target_length

    @property
    def query_start_distance(self) -> int:
        """Return bases before the aligned query interval."""
        return self.query_start

    @property
    def query_end_distance(self) -> int:
        """Return bases after the aligned query interval."""
        return self.query_length - self.query_end

    @property
    def target_start_distance(self) -> int:
        """Return bases before the aligned target interval in forward coordinates."""
        return self.target_start

    @property
    def target_end_distance(self) -> int:
        """Return bases after the aligned target interval in forward coordinates."""
        return self.target_length - self.target_end


@dataclass(frozen=True, slots=True)
class Relationship:
    """Explicit classification and diagnostics for one alignment."""

    relationship_type: RelationshipType
    query_id: str
    target_id: str
    orientation: Orientation
    identity: float
    aligned_length: int
    query_coverage: float
    target_coverage: float
    status: str
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class MergeComponent:
    """A planned graph component; it does not imply that merging is safe."""

    component_id: str
    sequence_ids: tuple[str, ...]
    relationship_ids: tuple[str, ...] = ()
    ambiguous: bool = False


@dataclass(frozen=True, slots=True)
class BaseEvidence:
    """Sample-specific evidence for a base on a source contig."""

    sample: str
    contig_id: str
    position: int
    allele_counts: dict[str, int]
    depth: int
    mean_base_quality: float | None = None
    mean_mapping_quality: float | None = None


@dataclass(frozen=True, slots=True)
class JoinEvidence:
    """Sample-specific evidence gathered for a proposed new junction."""

    sample: str
    left_contig_id: str
    right_contig_id: str
    spanning_reads: int | None = None
    testable: bool = False
    diagnostics: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ConsensusDecision:
    """A policy result kept separate from its underlying evidence."""

    position: int
    chosen_base: str | None
    policy: ConflictPolicy
    status: str
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RunConfig:
    """Normalised run configuration; percentages are stored as fractions."""

    identity: float = 0.98
    min_overlap: int = 1000
    min_containment: int = 500
    containment_coverage: float = 0.98
    end_tolerance: int = 50
    kmer_size: int = 21
    window_size: int = 10
    min_shared_minimisers: int = 5
    max_minimiser_frequency: int = 100
    threads: int = 1
    evidence: EvidenceMode = EvidenceMode.NONE
    conflict_policy: ConflictPolicy = ConflictPolicy.REJECT
    output_prefix: Path = Path("contigger")
    deterministic_seed: int | None = None
    emit_gfa: bool = False

    def __post_init__(self) -> None:
        if not 0.0 <= self.identity <= 1.0:
            raise ConfigurationError("identity must be between 0 and 100 percent")
        if not 0.0 <= self.containment_coverage <= 1.0:
            raise ConfigurationError("containment coverage must be between 0 and 100 percent")
        if self.min_overlap < 1 or self.min_containment < 1:
            raise ConfigurationError("minimum overlap and containment must be positive")
        if self.end_tolerance < 0:
            raise ConfigurationError("end tolerance cannot be negative")
        if self.kmer_size < 1 or self.window_size < 1:
            raise ConfigurationError("k-mer and minimiser window sizes must be positive")
        if self.min_shared_minimisers < 1 or self.max_minimiser_frequency < 1:
            raise ConfigurationError("minimiser thresholds must be positive")
        if self.threads < 1:
            raise ConfigurationError("thread count must be at least one")

    def as_dict(self) -> dict[str, Any]:
        """Return a deterministic JSON-compatible representation."""
        return {
            "conflict_policy": self.conflict_policy.value,
            "containment_coverage": self.containment_coverage,
            "deterministic_seed": self.deterministic_seed,
            "emit_gfa": self.emit_gfa,
            "end_tolerance": self.end_tolerance,
            "evidence": self.evidence.value,
            "identity": self.identity,
            "kmer_size": self.kmer_size,
            "max_minimiser_frequency": self.max_minimiser_frequency,
            "min_containment": self.min_containment,
            "min_overlap": self.min_overlap,
            "min_shared_minimisers": self.min_shared_minimisers,
            "output_prefix": str(self.output_prefix),
            "threads": self.threads,
            "window_size": self.window_size,
        }
