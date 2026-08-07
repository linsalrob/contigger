"""Stable typed data models shared across Contigger components."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
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


class GraphEdgeKind(StrEnum):
    """Structural role of a retained relationship in a sequence graph."""

    CONTAINMENT = "containment"
    OVERLAP = "overlap"
    AMBIGUOUS = "ambiguous"


class GraphDecisionStatus(StrEnum):
    """Conservative eligibility assigned by a graph decision policy."""

    ELIGIBLE = "eligible"
    DEFERRED = "deferred"


class JunctionSupportStatus(StrEnum):
    """Technology-scoped interpretation of targeted remapping evidence."""

    SUPPORTED = "supported"
    UNSUPPORTED = "unsupported"
    DEFERRED = "deferred"


class AlignmentType(StrEnum):
    """PAF alignment role reported by an aligner, when available."""

    PRIMARY = "P"
    SECONDARY = "S"
    SUPPLEMENTARY = "SUPPLEMENTARY"
    INVERSION_PRIMARY = "I"
    INVERSION_SECONDARY = "i"


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


class ContigEnd(StrEnum):
    """Named physical end of a source contig."""

    PREFIX = "prefix"
    SUFFIX = "suffix"


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
    query_positions: tuple[int, ...] = ()
    target_positions: tuple[int, ...] = ()
    supported_orientations: tuple[Orientation, ...] = ()
    terminal_topologies: tuple[str, ...] = ()
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class CatalogueSequence:
    """One canonical sequence retained after exact strand-aware deduplication."""

    identifier: str
    sequence: str
    length: int
    sha256: str
    representative_source_id: str

    def __post_init__(self) -> None:
        if self.length != len(self.sequence):
            raise InputValidationError("catalogue sequence length does not match its bases")


@dataclass(frozen=True, slots=True)
class CatalogueMember:
    """One source contig mapped recoverably onto a catalogue sequence."""

    catalogue_id: str
    source_id: str
    source_sample: str
    original_identifier: str
    orientation: Orientation
    representative: bool = False


@dataclass(frozen=True, slots=True)
class SequenceCatalogue:
    """Deterministically ordered canonical sequences and complete source membership."""

    sequences: tuple[CatalogueSequence, ...]
    members: tuple[CatalogueMember, ...]


@dataclass(frozen=True, slots=True)
class MinimiserObservation:
    """One canonical k-mer selected by a minimiser window at a known position."""

    sequence_id: str
    value: int
    position: int
    orientation: Orientation
    kmer: str


@dataclass(frozen=True, slots=True)
class AlignmentRequest:
    """Validated request to align exactly one candidate catalogue pair."""

    query: SequenceRecord
    target: SequenceRecord
    candidate: CandidatePair


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
    chaining_score: int | None = None
    secondary_chaining_score: int | None = None
    alignment_type: AlignmentType | None = None

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
        if self.mapping_quality is not None and not 0 <= self.mapping_quality <= 255:
            raise InputValidationError("mapping quality must be between 0 and 255")
        query_span = self.query_end - self.query_start
        target_span = self.target_end - self.target_start
        if self.matching_bases > min(query_span, target_span):
            raise InputValidationError("matching bases cannot exceed either aligned sequence span")
        if self.alignment_block_length < max(query_span, target_span):
            raise InputValidationError(
                "alignment block length cannot be shorter than an aligned sequence span"
            )
        if self.alignment_block_length > query_span + target_span:
            raise InputValidationError(
                "alignment block length cannot exceed the sum of aligned sequence spans"
            )

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
class RejectedAlignment:
    """A rejected alignment retained with its single-hit diagnostic."""

    hit: AlignmentHit
    relationship: Relationship


@dataclass(frozen=True, slots=True)
class PairRelationship:
    """Conservative decision over every distinct hit for one ordered pair."""

    relationship: Relationship
    representative_hit: AlignmentHit | None
    accepted_hits: tuple[AlignmentHit, ...]
    rejected_alignments: tuple[RejectedAlignment, ...]
    ambiguity_reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class GraphNode:
    """One stable sequence identity in a relationship graph."""

    sequence_id: str


@dataclass(frozen=True, slots=True)
class GraphEdge:
    """One canonical pairwise graph edge retaining alignment diagnostics."""

    edge_id: str
    kind: GraphEdgeKind
    relationship_type: RelationshipType
    query_id: str
    target_id: str
    orientation: Orientation
    query_start: int | None
    query_end: int | None
    target_start: int | None
    target_end: int | None
    identity: float
    aligned_length: int
    accepted_hit_count: int
    rejected_hit_count: int
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class MergeComponent:
    """A planned graph component; it does not imply that merging is safe."""

    component_id: str
    sequence_ids: tuple[str, ...]
    relationship_ids: tuple[str, ...] = ()
    ambiguous: bool = False
    ambiguity_reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RelationshipGraph:
    """Deterministic unsimplified graph with structurally separated edge classes."""

    nodes: tuple[GraphNode, ...]
    containment_edges: tuple[GraphEdge, ...]
    overlap_edges: tuple[GraphEdge, ...]
    ambiguous_edges: tuple[GraphEdge, ...]
    components: tuple[MergeComponent, ...]


@dataclass(frozen=True, slots=True)
class ContainmentDecision:
    """Policy decision for one containment edge without removing a sequence."""

    edge_id: str
    contained_sequence_id: str
    container_sequence_id: str
    status: GraphDecisionStatus
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class OverlapComponentDecision:
    """Eligibility of one overlap component for later path planning."""

    component_id: str
    sequence_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    status: GraphDecisionStatus
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class GraphDecisionPlan:
    """Deterministic policy output that authorizes no sequence modification."""

    containment_decisions: tuple[ContainmentDecision, ...]
    overlap_decisions: tuple[OverlapComponentDecision, ...]


@dataclass(frozen=True, slots=True)
class PlannedSourceMember:
    """One source contig retained in a path node with path-relative strand."""

    source_id: str
    source_sample: str
    original_identifier: str
    orientation: Orientation


@dataclass(frozen=True, slots=True)
class PlannedPathNode:
    """One catalogue sequence placed in an oriented linear path."""

    sequence_id: str
    orientation: Orientation
    source_members: tuple[PlannedSourceMember, ...]


@dataclass(frozen=True, slots=True)
class LinearPathPlan:
    """Canonical metadata-only plan for an eligible linear overlap component."""

    path_id: str
    component_id: str
    nodes: tuple[PlannedPathNode, ...]
    edge_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PathPlanningResult:
    """Deterministic path plans and every overlap component left deferred."""

    paths: tuple[LinearPathPlan, ...]
    deferred_component_ids: tuple[str, ...]


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
class JunctionRemappingRequest:
    """Request with a zero-based junction boundary inside one provisional reference."""

    sample: str
    left_contig_id: str
    left_end: ContigEnd
    right_contig_id: str
    right_end: ContigEnd
    provisional_reference: SequenceRecord
    junction_position: int
    extraction_distance: int = 1000
    minimum_spanning_flank: int = 20
    minimum_mapping_quality: int = 0

    def __post_init__(self) -> None:
        if not 0 < self.junction_position < self.provisional_reference.length:
            raise InputValidationError("junction position must be inside the provisional reference")
        if self.extraction_distance < 1:
            raise InputValidationError("read extraction distance must be positive")
        if self.minimum_spanning_flank < 1:
            raise InputValidationError("minimum spanning flank must be positive")
        if not 0 <= self.minimum_mapping_quality <= 255:
            raise InputValidationError("minimum mapping quality must be between 0 and 255")
        if self.junction_position < self.minimum_spanning_flank or (
            self.provisional_reference.length - self.junction_position < self.minimum_spanning_flank
        ):
            raise InputValidationError(
                "provisional reference is too short for the requested spanning flank"
            )


@dataclass(frozen=True, slots=True)
class TargetedJunctionEvidence:
    """Deterministic remapping observations that do not themselves authorize a merge."""

    sample: str
    technology: str
    remapping_preset: str
    left_contig_id: str
    right_contig_id: str
    provisional_reference_id: str
    provisional_reference_length: int
    provisional_reference_sha256: str
    junction_position: int
    selected_read_names: tuple[str, ...]
    remapped_read_names: tuple[str, ...]
    spanning_read_names: tuple[str, ...]
    minimum_spanning_flank: int
    samtools_version: str
    minimap2_version: str
    commands: tuple[tuple[str, ...], ...]
    diagnostics: tuple[str, ...] = ()
    minimum_mapping_quality: int = 0

    def __post_init__(self) -> None:
        if not self.technology or not self.remapping_preset:
            raise InputValidationError("junction evidence requires technology and remapping preset")
        if self.provisional_reference_length < 1:
            raise InputValidationError("provisional reference length must be positive")
        if len(self.provisional_reference_sha256) != 64 or any(
            symbol not in "0123456789abcdef" for symbol in self.provisional_reference_sha256
        ):
            raise InputValidationError(
                "provisional reference SHA-256 must be lowercase hexadecimal"
            )
        if not 0 < self.junction_position < self.provisional_reference_length:
            raise InputValidationError("junction position must be inside the provisional reference")
        if self.minimum_spanning_flank < 1 or (
            self.minimum_spanning_flank > self.junction_position
            or self.minimum_spanning_flank
            > self.provisional_reference_length - self.junction_position
        ):
            raise InputValidationError(
                "minimum spanning flank must fit on both sides of the junction"
            )
        if not 0 <= self.minimum_mapping_quality <= 255:
            raise InputValidationError("minimum mapping quality must be between 0 and 255")
        for label, names in (
            ("selected", self.selected_read_names),
            ("remapped", self.remapped_read_names),
            ("spanning", self.spanning_read_names),
        ):
            if len(names) != len(set(names)):
                raise InputValidationError(f"{label} junction read names must be unique")
        if not set(self.remapped_read_names) <= set(self.selected_read_names):
            raise InputValidationError("remapped reads must be a subset of selected reads")
        if not set(self.spanning_read_names) <= set(self.remapped_read_names):
            raise InputValidationError("spanning reads must be a subset of remapped reads")

    @property
    def spanning_reads(self) -> int:
        """Return the number of distinct reads spanning the provisional junction."""
        return len(self.spanning_read_names)


@dataclass(frozen=True, slots=True)
class JunctionPolicyReview:
    """External review provenance required before a junction policy is reviewed."""

    truth_dataset_sha256: str
    candidate_baseline_sha256: str
    reviewer: str
    reviewed_at: str
    decision: str
    technology: str
    remapping_preset: str
    minimum_spanning_reads: int
    minimum_spanning_fraction: float
    minimum_spanning_flank: int
    minimum_mapping_quality: int
    notes: str = ""

    def __post_init__(self) -> None:
        for label, value in (
            ("truth dataset", self.truth_dataset_sha256),
            ("candidate baseline", self.candidate_baseline_sha256),
        ):
            if len(value) != 64 or any(symbol not in "0123456789abcdef" for symbol in value):
                raise InputValidationError(f"{label} SHA-256 must be lowercase hexadecimal")
        if not self.reviewer.strip():
            raise InputValidationError("junction policy review requires reviewer and timestamp")
        try:
            parsed_timestamp = datetime.fromisoformat(self.reviewed_at.replace("Z", "+00:00"))
        except ValueError as error:
            raise InputValidationError(
                "junction policy review timestamp must be RFC 3339"
            ) from error
        if "T" not in self.reviewed_at or parsed_timestamp.tzinfo is None:
            raise InputValidationError("junction policy review timestamp must be RFC 3339")
        if self.decision not in {"approved", "rejected"}:
            raise InputValidationError(
                "junction policy review decision must be approved or rejected"
            )
        if not self.technology or not self.remapping_preset:
            raise InputValidationError("junction policy review requires technology and preset")
        if self.minimum_spanning_reads < 1 or self.minimum_spanning_flank < 1:
            raise InputValidationError("junction policy review thresholds must be positive")
        if not 0.0 <= self.minimum_spanning_fraction <= 1.0:
            raise InputValidationError(
                "junction policy review fraction must be between zero and one"
            )
        if not 0 <= self.minimum_mapping_quality <= 255:
            raise InputValidationError(
                "junction policy review mapping quality must be between 0 and 255"
            )


@dataclass(frozen=True, slots=True)
class JunctionTruthSetMetadata:
    """Auditable scope and counts for one technology-specific junction truth set."""

    dataset_id: str
    dataset_version: str
    technology: str
    remapping_preset: str
    source_description: str
    truth_sha256: str
    case_count: int
    true_case_count: int
    artificial_case_count: int
    false_support_baseline_established: bool
    reviewed: bool = False

    def __post_init__(self) -> None:
        if not all(
            value.strip()
            for value in (
                self.dataset_id,
                self.dataset_version,
                self.technology,
                self.remapping_preset,
                self.source_description,
            )
        ):
            raise InputValidationError("junction truth metadata identifiers cannot be blank")
        if len(self.truth_sha256) != 64 or any(
            symbol not in "0123456789abcdef" for symbol in self.truth_sha256
        ):
            raise InputValidationError("junction truth SHA-256 must be lowercase hexadecimal")
        counts = (self.case_count, self.true_case_count, self.artificial_case_count)
        if any(type(value) is not int or value < 0 for value in counts):
            raise InputValidationError(
                "junction truth metadata counts must be non-negative integers"
            )
        if self.true_case_count + self.artificial_case_count != self.case_count:
            raise InputValidationError("junction truth metadata case counts do not balance")
        if (
            type(self.false_support_baseline_established) is not bool
            or type(self.reviewed) is not bool
        ):
            raise InputValidationError("junction truth metadata flags must be Boolean")
        if self.reviewed and not self.false_support_baseline_established:
            raise InputValidationError(
                "reviewed junction truth requires an established false-support baseline"
            )


@dataclass(frozen=True, slots=True)
class JunctionSupportPolicy:
    """Explicit technology-specific criteria requiring external benchmark review."""

    technology: str
    remapping_preset: str
    minimum_spanning_reads: int
    minimum_spanning_fraction: float
    minimum_spanning_flank: int
    reviewed: bool = False
    minimum_mapping_quality: int = 0
    review: JunctionPolicyReview | None = None
    truth_metadata: JunctionTruthSetMetadata | None = None
    candidate_baseline_sha256: str | None = None

    def __post_init__(self) -> None:
        if not self.technology or not self.remapping_preset:
            raise InputValidationError(
                "junction support policy requires a technology and remapping preset"
            )
        if self.minimum_spanning_reads < 1 or self.minimum_spanning_flank < 1:
            raise InputValidationError("junction support counts and flank must be positive")
        if not 0.0 <= self.minimum_spanning_fraction <= 1.0:
            raise InputValidationError("minimum spanning fraction must be between zero and one")
        if not 0 <= self.minimum_mapping_quality <= 255:
            raise InputValidationError("minimum mapping quality must be between 0 and 255")
        if self.reviewed and (self.review is None or self.review.decision != "approved"):
            raise InputValidationError(
                "reviewed junction policy requires an approved review artifact"
            )
        if self.reviewed and (
            self.truth_metadata is None or self.candidate_baseline_sha256 is None
        ):
            raise InputValidationError(
                "reviewed junction policy requires paired truth metadata and candidate baseline"
            )
        if self.reviewed and self.review is not None:
            configuration = (
                self.technology,
                self.remapping_preset,
                self.minimum_spanning_reads,
                self.minimum_spanning_fraction,
                self.minimum_spanning_flank,
                self.minimum_mapping_quality,
            )
            reviewed_configuration = (
                self.review.technology,
                self.review.remapping_preset,
                self.review.minimum_spanning_reads,
                self.review.minimum_spanning_fraction,
                self.review.minimum_spanning_flank,
                self.review.minimum_mapping_quality,
            )
            if configuration != reviewed_configuration:
                raise InputValidationError(
                    "reviewed junction policy does not match its approved review artifact"
                )
        if self.reviewed and self.review is not None and self.truth_metadata is not None:
            if (
                not self.truth_metadata.reviewed
                or not self.truth_metadata.false_support_baseline_established
            ):
                raise InputValidationError(
                    "reviewed junction policy requires reviewed truth metadata and "
                    "false-support baseline"
                )
            if self.review.truth_dataset_sha256 != self.truth_metadata.truth_sha256:
                raise InputValidationError(
                    "reviewed junction policy truth digest does not match metadata"
                )
            if (
                self.review.technology != self.truth_metadata.technology
                or self.review.remapping_preset != self.truth_metadata.remapping_preset
            ):
                raise InputValidationError(
                    "reviewed junction policy review technology or preset does not "
                    "match truth metadata"
                )
            if self.review.candidate_baseline_sha256 != self.candidate_baseline_sha256:
                raise InputValidationError(
                    "reviewed junction policy candidate baseline digest does not "
                    "match supplied baseline"
                )


@dataclass(frozen=True, slots=True)
class JunctionSupportDecision:
    """Evidence interpretation kept separate from graph eligibility and merging."""

    status: JunctionSupportStatus
    technology: str
    spanning_reads: int
    remapped_reads: int
    spanning_fraction: float
    reasons: tuple[str, ...]


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
            raise ConfigurationError("identity must be a fraction between 0.0 and 1.0")
        if not 0.0 <= self.containment_coverage <= 1.0:
            raise ConfigurationError("containment_coverage must be a fraction between 0.0 and 1.0")
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
