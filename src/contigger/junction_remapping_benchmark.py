"""Deterministic checked-in FASTQ remapping benchmark for proposed junctions."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
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


def evaluate_junction_remapping_benchmark(
    dataset: Path,
    *,
    minimap2: str | Path = "minimap2",
    preset: str = "map-ont",
    threads: int = 1,
    minimum_spanning_flank: int = 20,
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
    for expected in truth:
        reports = []
        details = []
        for sample in samples:
            fastq = dataset / "reads" / f"{sample.sample}.targeted.fastq.gz"
            if not fastq.is_file():
                raise InputValidationError(f"junction benchmark FASTQ does not exist: {fastq}")
            request = requests[expected.case_id]
            request = JunctionRemappingRequest(
                sample.sample,
                request.left_contig_id,
                request.left_end,
                request.right_contig_id,
                request.right_end,
                request.provisional_reference,
                request.junction_position,
                minimum_spanning_flank=minimum_spanning_flank,
            )
            report = remapper.evaluate(
                request,
                sample=sample.sample,
                technology=sample.technology or "unknown",
                fastq=fastq,
            )
            reports.append(report)
            details.append(
                {
                    "sample": report.sample,
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
            f"false spanning-support sample cases: {summary.false_support_sample_cases}",
            f"missed spanning-support sample cases: {summary.missed_support_sample_cases}",
            "false-support baseline established: "
            + ("yes" if summary.false_support_baseline_established else "no"),
            "testable artificial sample cases: "
            f"{summary.testable_artificial_sample_cases}/{summary.artificial_sample_cases}",
            f"dataset version: {report.dataset_version}",
            f"preset: {report.preset}",
            f"minimum spanning flank: {report.minimum_spanning_flank}",
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


def _names_sha256(names: tuple[str, ...]) -> str:
    return hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()
