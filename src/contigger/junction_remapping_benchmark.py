"""Deterministic checked-in FASTQ remapping benchmark for proposed junctions."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, TextIO

from contigger import __version__
from contigger.benchmark import load_truth
from contigger.evidence.junctions import FastqJunctionRemapper
from contigger.exceptions import InputValidationError
from contigger.fasta import read_fasta
from contigger.junction_benchmark import load_junction_truth
from contigger.manifest import parse_manifest
from contigger.models import (
    ContigEnd,
    JunctionRemappingRequest,
    Orientation,
    RelationshipType,
    SequenceRecord,
    TargetedJunctionEvidence,
)


@dataclass(frozen=True, slots=True)
class JunctionRemappingBenchmarkReport:
    """Reproducible metadata and scored sample-specific remapping observations."""

    dataset_version: str
    dataset_truth_sha256: str
    preset: str
    minimum_spanning_flank: int
    minimum_mapping_quality: int
    minimap2_version: str
    contigger_version: str
    summary: SampleScopedJunctionSummary
    cases: tuple[SampleScopedJunctionCase, ...]
    aggregate_case_status: dict[str, tuple[str, ...]]
    observations: dict[str, tuple[dict[str, Any], ...]]


@dataclass(frozen=True, slots=True)
class SampleScopedJunctionCase:
    """One truth case evaluated in one sample without pooling read counts."""

    case_id: str
    sample: str
    observation_kind: str
    junction_is_true: bool
    remapped_reads: int
    spanning_reads: int
    testable: bool
    false_support: bool
    missed_support: bool


@dataclass(frozen=True, slots=True)
class SampleScopedJunctionSummary:
    """Counts over independent sample-case observations."""

    expected_junctions: int
    expected_sample_cases: int
    true_sample_cases: int
    artificial_sample_cases: int
    testable_artificial_sample_cases: int
    false_support_sample_cases: int
    missed_support_sample_cases: int
    false_support_baseline_established: bool


@dataclass(frozen=True, slots=True)
class JunctionPolicyCandidateBenchmark:
    """Score one unreviewed support-threshold candidate against observations."""

    minimum_spanning_reads: int
    minimum_spanning_fraction: float
    false_support_cases: int
    missed_true_cases: int
    supported_true_cases: int
    testable_negative_cases: int
    reviewed: bool = False


@dataclass(frozen=True, slots=True)
class JunctionPolicyCandidateBaselineArtifact:
    """Strictly parsed candidate benchmark provenance, without policy authorization."""

    candidates: tuple[JunctionPolicyCandidateBenchmark, ...]
    negative_controls: int
    reviewed: bool
    result: str


def load_junction_policy_candidate_baseline(
    path: Path, *, expected_sha256: str | None = None
) -> JunctionPolicyCandidateBaselineArtifact:
    """Load a candidate baseline and optionally verify its exact file digest."""
    try:
        raw = path.read_bytes()
        payload = json.loads(raw.decode("utf-8"), object_pairs_hook=_reject_duplicate_keys)
    except (OSError, UnicodeError, ValueError) as error:
        raise InputValidationError(
            f"cannot read junction policy candidate baseline {path}: {error}"
        ) from error
    if expected_sha256 is not None and hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise InputValidationError(
            f"{path}: candidate baseline digest does not match expected value"
        )
    if not isinstance(payload, dict):
        raise InputValidationError(f"{path}: candidate baseline must be a JSON object")
    required = {"candidates", "negative_controls", "reviewed", "result"}
    missing = sorted(required - set(payload))
    unknown = sorted(set(payload) - required)
    if missing:
        raise InputValidationError(f"{path}: candidate baseline is missing {missing[0]}")
    if unknown:
        raise InputValidationError(f"{path}: candidate baseline has unknown field {unknown[0]!r}")
    if type(payload["negative_controls"]) is not int or payload["negative_controls"] < 0:
        raise InputValidationError(f"{path}: negative_controls must be a non-negative integer")
    if type(payload["reviewed"]) is not bool or not isinstance(payload["result"], str):
        raise InputValidationError(f"{path}: invalid candidate baseline metadata")
    if payload["reviewed"]:
        raise InputValidationError(
            f"{path}: candidate baseline reviewed status requires an external policy review"
        )
    candidates = payload["candidates"]
    if not isinstance(candidates, list) or not candidates:
        raise InputValidationError(f"{path}: candidates must be a non-empty array")
    candidate_fields = {
        "minimum_spanning_fraction",
        "minimum_spanning_reads",
        "false_support_cases",
        "missed_true_cases",
        "supported_true_cases",
    }
    parsed: list[JunctionPolicyCandidateBenchmark] = []
    seen: set[tuple[int, float]] = set()
    true_case_total: int | None = None
    for index, candidate in enumerate(candidates):
        if not isinstance(candidate, dict) or set(candidate) != candidate_fields:
            raise InputValidationError(f"{path}: candidate {index} has invalid fields")
        if (
            type(candidate["minimum_spanning_reads"]) is not int
            or candidate["minimum_spanning_reads"] < 1
        ):
            raise InputValidationError(
                f"{path}: candidate {index} has invalid minimum_spanning_reads"
            )
        if (
            type(candidate["minimum_spanning_fraction"]) not in {int, float}
            or not 0.0 <= candidate["minimum_spanning_fraction"] <= 1.0
        ):
            raise InputValidationError(
                f"{path}: candidate {index} has invalid minimum_spanning_fraction"
            )
        counts = tuple(
            candidate[field]
            for field in ("false_support_cases", "missed_true_cases", "supported_true_cases")
        )
        if any(type(value) is not int or value < 0 for value in counts):
            raise InputValidationError(f"{path}: candidate {index} has invalid result counts")
        if counts[0] > payload["negative_controls"]:
            raise InputValidationError(
                f"{path}: candidate {index} false support exceeds negative controls"
            )
        candidate_true_total = counts[1] + counts[2]
        if true_case_total is None:
            true_case_total = candidate_true_total
        elif candidate_true_total != true_case_total:
            raise InputValidationError(
                f"{path}: candidate {index} true-case totals are inconsistent"
            )
        key = (candidate["minimum_spanning_reads"], float(candidate["minimum_spanning_fraction"]))
        if key in seen:
            raise InputValidationError(f"{path}: duplicate candidate threshold at index {index}")
        seen.add(key)
        parsed.append(
            JunctionPolicyCandidateBenchmark(
                minimum_spanning_reads=key[0],
                minimum_spanning_fraction=key[1],
                false_support_cases=counts[0],
                missed_true_cases=counts[1],
                supported_true_cases=counts[2],
                testable_negative_cases=payload["negative_controls"],
                reviewed=payload["reviewed"],
            )
        )
    return JunctionPolicyCandidateBaselineArtifact(
        tuple(parsed), payload["negative_controls"], payload["reviewed"], payload["result"]
    )


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON object members instead of silently choosing one."""
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field {key!r}")
        result[key] = value
    return result


def benchmark_junction_policy_candidates(
    report: JunctionRemappingBenchmarkReport,
    candidates: tuple[tuple[int, float], ...],
) -> tuple[JunctionPolicyCandidateBenchmark, ...]:
    """Evaluate threshold candidates without authorizing evidence or graph edges."""
    results: list[JunctionPolicyCandidateBenchmark] = []
    for minimum_reads, minimum_fraction in candidates:
        if minimum_reads < 1:
            raise InputValidationError("candidate minimum spanning reads must be positive")
        if not 0.0 <= minimum_fraction <= 1.0:
            raise InputValidationError("candidate spanning fraction must be between zero and one")
        false_support = 0
        missed_true = 0
        supported_true = 0
        testable_negative = 0
        for case in report.cases:
            fraction = case.spanning_reads / case.remapped_reads if case.remapped_reads else 0.0
            supported = case.spanning_reads >= minimum_reads and fraction >= minimum_fraction
            if not case.junction_is_true:
                testable_negative += case.testable
                false_support += supported
            elif supported:
                supported_true += 1
            else:
                missed_true += 1
        results.append(
            JunctionPolicyCandidateBenchmark(
                minimum_reads,
                minimum_fraction,
                false_support,
                missed_true,
                supported_true,
                testable_negative,
            )
        )
    return tuple(results)


def evaluate_junction_remapping_benchmark(
    dataset: Path,
    *,
    minimap2: str | Path = "minimap2",
    preset: str = "map-ont",
    threads: int = 1,
    minimum_spanning_flank: int = 20,
    minimum_mapping_quality: int = 0,
) -> JunctionRemappingBenchmarkReport:
    """Remap every checked-in sample FASTQ to every explicit benchmark junction."""
    required = (
        dataset / "VERSION",
        dataset / "manifest.tsv",
        dataset / "expected" / "expected_junctions.tsv",
        dataset / "expected" / "expected_relationships.tsv",
        dataset / "expected" / "expected_merged_sequences.fasta.gz",
    )
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise InputValidationError(f"junction benchmark dataset is missing {missing[0]}")
    version = (dataset / "VERSION").read_text(encoding="utf-8").strip()
    if not version:
        raise InputValidationError("junction benchmark dataset VERSION is empty")
    samples = parse_manifest(dataset / "manifest.tsv").samples
    truth_path = dataset / "expected" / "expected_junctions.tsv"
    truth = load_junction_truth(truth_path)
    requests = _build_requests(dataset, minimum_spanning_flank)
    remapper = FastqJunctionRemapper(minimap2=minimap2, preset=preset, threads=threads)
    evidence_by_case: dict[str, tuple[TargetedJunctionEvidence, ...]] = {}
    serialized: dict[str, tuple[dict[str, Any], ...]] = {}
    with TemporaryDirectory(prefix="contigger-negative-controls-") as directory:
        controls = _write_artificial_negative_controls(dataset, Path(directory))
        for expected in truth:
            reports = []
            details = []
            inputs = (
                tuple(
                    (
                        sample.sample,
                        sample.technology or "unknown",
                        dataset / "reads" / f"{sample.sample}.targeted.fastq.gz",
                        "sample-targeted-reads",
                    )
                    for sample in samples
                )
                if expected.junction_is_true
                else (
                    (
                        "synthetic-negative-control",
                        "ont",
                        controls[expected.case_id],
                        "synthetic-source-end-control",
                    ),
                )
            )
            for sample_id, technology, fastq, observation_kind in inputs:
                if not fastq.is_file():
                    raise InputValidationError(f"junction benchmark FASTQ does not exist: {fastq}")
                request = requests[expected.case_id]
                request = JunctionRemappingRequest(
                    sample_id,
                    request.left_contig_id,
                    request.left_end,
                    request.right_contig_id,
                    request.right_end,
                    request.provisional_reference,
                    request.junction_position,
                    minimum_spanning_flank=minimum_spanning_flank,
                    minimum_mapping_quality=minimum_mapping_quality,
                )
                report = remapper.evaluate(
                    request,
                    sample=sample_id,
                    technology=technology,
                    fastq=fastq,
                )
                reports.append(report)
                details.append(
                    {
                        "sample": report.sample,
                        "observation_kind": observation_kind,
                        "selected_reads": len(report.selected_read_names),
                        "remapped_reads": len(report.remapped_read_names),
                        "spanning_reads": report.spanning_reads,
                        "spanning_read_names_sha256": _names_sha256(report.spanning_read_names),
                        "provisional_reference_sha256": report.provisional_reference_sha256,
                    }
                )
            evidence_by_case[expected.case_id] = tuple(reports)
            serialized[expected.case_id] = tuple(details)
    cases = tuple(
        SampleScopedJunctionCase(
            expected.case_id,
            report.sample,
            "sample-targeted-reads"
            if expected.junction_is_true
            else "synthetic-source-end-control",
            expected.junction_is_true,
            len(report.remapped_read_names),
            report.spanning_reads,
            bool(report.remapped_read_names),
            not expected.junction_is_true and report.spanning_reads > 0,
            expected.junction_is_true and report.spanning_reads == 0,
        )
        for expected in truth
        for report in evidence_by_case[expected.case_id]
    )
    artificial = tuple(case for case in cases if not case.junction_is_true)
    testable_artificial = tuple(case for case in artificial if case.testable)
    summary = SampleScopedJunctionSummary(
        expected_junctions=len(truth),
        expected_sample_cases=len(cases),
        true_sample_cases=sum(case.junction_is_true for case in cases),
        artificial_sample_cases=len(artificial),
        testable_artificial_sample_cases=len(testable_artificial),
        false_support_sample_cases=sum(case.false_support for case in cases),
        missed_support_sample_cases=sum(case.missed_support for case in cases),
        false_support_baseline_established=len(testable_artificial) == len(artificial),
    )
    aggregate_status = {
        "detected_true_cases": tuple(
            expected.case_id
            for expected in truth
            if expected.junction_is_true
            and any(report.spanning_reads > 0 for report in evidence_by_case[expected.case_id])
        ),
        "missed_true_cases": tuple(
            expected.case_id
            for expected in truth
            if expected.junction_is_true
            and not any(report.spanning_reads > 0 for report in evidence_by_case[expected.case_id])
        ),
        "artificial_cases_with_support": tuple(
            expected.case_id
            for expected in truth
            if not expected.junction_is_true
            and any(report.spanning_reads > 0 for report in evidence_by_case[expected.case_id])
        ),
        "untestable_artificial_cases": tuple(
            expected.case_id
            for expected in truth
            if not expected.junction_is_true
            and not any(report.remapped_read_names for report in evidence_by_case[expected.case_id])
        ),
    }
    minimap_version = next(
        report.minimap2_version for reports in evidence_by_case.values() for report in reports
    )
    return JunctionRemappingBenchmarkReport(
        version,
        hashlib.sha256(truth_path.read_bytes()).hexdigest(),
        preset,
        minimum_spanning_flank,
        minimum_mapping_quality,
        minimap_version,
        __version__,
        summary,
        cases,
        aggregate_status,
        serialized,
    )


def write_junction_remapping_json(report: JunctionRemappingBenchmarkReport, output: TextIO) -> None:
    """Write deterministic JSON without volatile temporary command paths."""
    json.dump(asdict(report), output, indent=2, sort_keys=True)
    output.write("\n")


def format_junction_remapping_summary(report: JunctionRemappingBenchmarkReport) -> str:
    """Format conservative false-support-first benchmark results."""
    summary = report.summary
    return "\n".join(
        (
            f"false spanning-support observation cases: {summary.false_support_sample_cases}",
            f"missed spanning-support observation cases: {summary.missed_support_sample_cases}",
            "false-support baseline established: "
            + ("yes" if summary.false_support_baseline_established else "no"),
            "testable artificial controls: "
            f"{summary.testable_artificial_sample_cases}/{summary.artificial_sample_cases}",
            f"dataset version: {report.dataset_version}",
            f"preset: {report.preset}",
            f"minimum spanning flank: {report.minimum_spanning_flank}",
            f"minimum mapping quality: {report.minimum_mapping_quality}",
            f"expected junctions: {summary.expected_junctions}",
            "result: remapping evidence only; no merge was authorized",
        )
    )


def _build_requests(dataset: Path, flank: int) -> dict[str, JunctionRemappingRequest]:
    samples = parse_manifest(dataset / "manifest.tsv").samples
    sequences: dict[str, SequenceRecord] = {}
    for sample in samples:
        for record in read_fasta(sample.contigs):
            if record.identifier in sequences:
                raise InputValidationError(
                    f"duplicate source contig identifier {record.identifier!r} in benchmark"
                )
            sequences[record.identifier] = record
    relationships = load_truth(dataset / "expected" / "expected_relationships.tsv")
    expected_sequences = {
        item.identifier: item.sequence
        for item in read_fasta(dataset / "expected" / "expected_merged_sequences.fasta.gz")
    }
    requests: dict[str, JunctionRemappingRequest] = {}
    for truth in load_junction_truth(dataset / "expected" / "expected_junctions.tsv"):
        matches = [
            item
            for item in relationships
            if item.case_id == truth.case_id
            and item.query_id == truth.left_contig_id
            and item.target_id == truth.right_contig_id
        ]
        if len(matches) != 1:
            raise InputValidationError(
                f"junction case {truth.case_id!r} requires exactly one ordered relationship"
            )
        relationship = matches[0]
        if (
            relationship.relationship_type is not RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
            or relationship.orientation is not Orientation.FORWARD
        ):
            raise InputValidationError(
                f"junction case {truth.case_id!r} is not a forward suffix-to-prefix relationship"
            )
        try:
            left = sequences[truth.left_contig_id]
            right = sequences[truth.right_contig_id]
        except KeyError as error:
            raise InputValidationError(
                f"junction case {truth.case_id!r} references an unknown source contig"
            ) from error
        sequence = left.sequence + right.sequence[relationship.target_end :]
        if truth.junction_is_true and expected_sequences.get(truth.case_id) != sequence:
            raise InputValidationError(
                f"junction case {truth.case_id!r} does not match expected merged sequence"
            )
        if not truth.junction_is_true and truth.case_id in expected_sequences:
            raise InputValidationError(
                f"artificial junction case {truth.case_id!r} has an expected merged sequence"
            )
        reference = SequenceRecord(
            identifier=f"junction_{truth.case_id}",
            source_sample="benchmark",
            original_identifier=truth.case_id,
            description="benchmark-only provisional junction",
            sequence=sequence,
            length=len(sequence),
        )
        requests[truth.case_id] = JunctionRemappingRequest(
            "benchmark",
            truth.left_contig_id,
            ContigEnd.SUFFIX,
            truth.right_contig_id,
            ContigEnd.PREFIX,
            reference,
            left.length,
            minimum_spanning_flank=flank,
        )
    return requests


def _write_artificial_negative_controls(dataset: Path, directory: Path) -> dict[str, Path]:
    """Create exact, non-spanning source-end reads for each artificial adjacency."""
    sequences: dict[str, str] = {}
    for sample in parse_manifest(dataset / "manifest.tsv").samples:
        for record in read_fasta(sample.contigs):
            if record.identifier in sequences:
                raise InputValidationError(
                    f"duplicate source contig identifier {record.identifier!r} in benchmark"
                )
            sequences[record.identifier] = record.sequence
    relationships = {
        item.case_id: item
        for item in load_truth(dataset / "expected" / "expected_relationships.tsv")
        if item.relationship_type is RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX
        and item.orientation is Orientation.FORWARD
    }
    controls: dict[str, Path] = {}
    for expected in load_junction_truth(dataset / "expected" / "expected_junctions.tsv"):
        if expected.junction_is_true:
            continue
        try:
            relationship = relationships[expected.case_id]
            left = sequences[expected.left_contig_id]
            right = sequences[expected.right_contig_id]
        except KeyError as error:
            raise InputValidationError(
                f"artificial junction {expected.case_id!r} lacks source-end control inputs"
            ) from error
        read_length = min(1000, len(left), len(right) - relationship.target_end)
        if read_length < 100:
            raise InputValidationError(
                f"artificial junction {expected.case_id!r} cannot supply 100 bp controls"
            )
        records = (
            (f"{expected.case_id}:left-end", left[-read_length:]),
            (
                f"{expected.case_id}:right-after-overlap",
                right[relationship.target_end : relationship.target_end + read_length],
            ),
        )
        path = directory / f"{expected.case_id}.fastq"
        path.write_text(
            "".join(
                f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n" for name, sequence in records
            ),
            encoding="ascii",
        )
        controls[expected.case_id] = path
    return controls


def _names_sha256(names: tuple[str, ...]) -> str:
    return hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()
