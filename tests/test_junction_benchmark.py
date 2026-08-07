"""Checked-in junction truth and conservative support-policy tests."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from contigger.exceptions import InputValidationError
from contigger.junction_benchmark import (
    evaluate_junction_support,
    load_junction_policy_review,
    load_junction_truth,
    load_junction_truth_set_metadata,
    score_junction_observations,
    validate_junction_policy_review,
)
from contigger.models import (
    JunctionPolicyReview,
    JunctionSupportPolicy,
    JunctionSupportStatus,
    JunctionTruthSetMetadata,
    TargetedJunctionEvidence,
)

DATASET = Path(__file__).parents[1] / "test_data"
REVIEW = JunctionPolicyReview(
    "a" * 64,
    "b" * 64,
    "reviewer",
    "2026-08-07T00:00:00Z",
    "approved",
    "ont",
    "map-ont",
    3,
    0.3,
    20,
    0,
)
METADATA = JunctionTruthSetMetadata(
    "dataset",
    "1",
    "ont",
    "map-ont",
    "synthetic reviewed truth",
    "a" * 64,
    19,
    15,
    4,
    True,
    True,
)


def review_for(policy: JunctionSupportPolicy) -> JunctionPolicyReview:
    """Build approved test provenance matching one policy configuration."""
    return replace(
        REVIEW,
        minimum_spanning_reads=policy.minimum_spanning_reads,
        minimum_spanning_fraction=policy.minimum_spanning_fraction,
        minimum_spanning_flank=policy.minimum_spanning_flank,
        minimum_mapping_quality=policy.minimum_mapping_quality,
    )


def evidence(case: str, spanning: int, remapped: int = 10) -> TargetedJunctionEvidence:
    """Build deterministic evidence with distinct sample-scoped read identities."""
    return TargetedJunctionEvidence(
        sample="S01",
        technology="ont",
        remapping_preset="map-ont",
        left_contig_id=f"{case}_left",
        right_contig_id=f"{case}_right",
        provisional_reference_id=case,
        provisional_reference_length=100,
        provisional_reference_sha256="a" * 64,
        junction_position=50,
        selected_read_names=tuple(f"read-{i}" for i in range(remapped)),
        remapped_read_names=tuple(f"read-{i}" for i in range(remapped)),
        spanning_read_names=tuple(f"read-{i}" for i in range(spanning)),
        minimum_spanning_flank=20,
        minimum_mapping_quality=0,
        samtools_version="samtools test",
        minimap2_version="minimap2 test",
        commands=(),
    )


def test_checked_in_junction_truth_is_typed_complete_and_sorted() -> None:
    truth = load_junction_truth(DATASET / "expected" / "expected_junctions.tsv")
    assert len(truth) == 19
    assert sum(item.junction_is_true for item in truth) == 15
    assert sum(not item.junction_is_true for item in truth) == 4
    assert sum(item.selected_spanning_reads for item in truth) == 1500
    assert [item.case_id for item in truth] == sorted(item.case_id for item in truth)
    circular = next(item for item in truth if item.case_id == "circular_origin")
    assert circular.circular_wrap and circular.source_coordinate == 0


def test_checked_in_junction_truth_baseline_matches_source() -> None:
    path = DATASET / "expected" / "expected_junctions.tsv"
    baseline = json.loads(
        (
            Path(__file__).parents[1] / "benchmarks" / "pseudomonas_junction_truth_baseline.json"
        ).read_text(encoding="utf-8")
    )
    assert baseline["junction_cases"] == 19
    assert baseline["true_junctions"] == 15
    assert baseline["artificial_junctions"] == 4
    assert baseline["expected_junctions_checksum"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_truth_errors_include_physical_line(tmp_path: Path) -> None:
    source = DATASET / "expected" / "expected_junctions.tsv"
    lines = source.read_text(encoding="utf-8").splitlines()
    path = tmp_path / "junctions.tsv"
    path.write_text("\n".join((lines[0], lines[1], lines[1])) + "\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match=r"junctions\.tsv:3: duplicate junction"):
        load_junction_truth(path)


def test_unreviewed_policy_is_always_deferred() -> None:
    policy = JunctionSupportPolicy("ont", "map-ont", 3, 0.5, 20)
    decision = evaluate_junction_support((evidence("case", 10),), policy)
    assert decision.status is JunctionSupportStatus.DEFERRED
    assert "not been reviewed" in decision.reasons[0]
    absent = evaluate_junction_support(
        (),
        replace(
            policy,
            reviewed=True,
            review=review_for(policy),
            truth_metadata=METADATA,
            candidate_baseline_sha256="b" * 64,
        ),
    )
    assert absent.status is JunctionSupportStatus.DEFERRED
    assert "no targeted-remapping" in absent.reasons[0]


def test_reviewed_policy_applies_inclusive_thresholds_and_flank_guard() -> None:
    policy = JunctionSupportPolicy(
        "ont",
        "map-ont",
        3,
        0.3,
        20,
        reviewed=True,
        review=REVIEW,
        truth_metadata=METADATA,
        candidate_baseline_sha256="b" * 64,
    )
    assert (
        evaluate_junction_support((evidence("case", 3),), policy).status
        is JunctionSupportStatus.SUPPORTED
    )
    assert (
        evaluate_junction_support((evidence("case", 2),), policy).status
        is JunctionSupportStatus.UNSUPPORTED
    )
    short_flank = replace(evidence("case", 10), minimum_spanning_flank=19)
    assert (
        evaluate_junction_support((short_flank,), policy).status is JunctionSupportStatus.DEFERRED
    )
    filtered = replace(evidence("case", 10), minimum_mapping_quality=20)
    mapping_quality_mismatch = evaluate_junction_support((filtered,), policy)
    assert mapping_quality_mismatch.status is JunctionSupportStatus.DEFERRED
    assert "mapping-quality" in mapping_quality_mismatch.reasons[0]
    matching_policy = replace(
        policy,
        minimum_mapping_quality=20,
        review=replace(REVIEW, minimum_mapping_quality=20),
    )
    assert (
        evaluate_junction_support((filtered,), matching_policy).status
        is JunctionSupportStatus.SUPPORTED
    )


def test_observational_scoring_separates_false_and_missed_support() -> None:
    truth = load_junction_truth(DATASET / "expected" / "expected_junctions.tsv")
    observations = {
        "circular_origin": (evidence("circular_origin", 1),),
        "end_tolerance_51": (evidence("end_tolerance_51", 1),),
    }
    report = score_junction_observations(truth, observations)
    assert report.summary.expected_cases == 19
    assert report.summary.observed_cases == 2
    assert report.summary.false_support_cases == 1
    assert report.summary.missed_support_cases == 14
    assert report.summary.correct_detections == 4
    assert report.cases == tuple(sorted(report.cases, key=lambda item: item.case_id))


def test_observations_reject_unknown_cases() -> None:
    truth = load_junction_truth(DATASET / "expected" / "expected_junctions.tsv")
    with pytest.raises(InputValidationError, match="unknown case"):
        score_junction_observations(truth, {"invented": (evidence("invented", 1),)})


def test_evidence_groups_reject_duplicate_samples_and_inconsistent_references() -> None:
    policy = JunctionSupportPolicy(
        "ont",
        "map-ont",
        1,
        0.0,
        20,
        reviewed=True,
        review=review_for(JunctionSupportPolicy("ont", "map-ont", 1, 0.0, 20)),
        truth_metadata=METADATA,
        candidate_baseline_sha256="b" * 64,
    )
    first = evidence("case", 1)
    with pytest.raises(InputValidationError, match="duplicate sample"):
        evaluate_junction_support((first, first), policy)
    second = replace(first, sample="S02", provisional_reference_sha256="b" * 64)
    with pytest.raises(InputValidationError, match="identical provisional junction"):
        evaluate_junction_support((first, second), policy)


def test_policy_never_pools_samples_and_requires_matching_configuration() -> None:
    policy = JunctionSupportPolicy(
        "ont",
        "map-ont",
        2,
        0.0,
        20,
        reviewed=True,
        review=review_for(JunctionSupportPolicy("ont", "map-ont", 2, 0.0, 20)),
        truth_metadata=METADATA,
        candidate_baseline_sha256="b" * 64,
    )
    first = replace(evidence("case", 1), remapped_read_names=("read-0",))
    second = replace(
        first,
        sample="S02",
        selected_read_names=("other",),
        remapped_read_names=("other",),
        spanning_read_names=("other",),
    )
    pooled = evaluate_junction_support((first, second), policy)
    assert pooled.status is JunctionSupportStatus.DEFERRED
    assert "cross-sample" in pooled.reasons[0]
    mismatch = evaluate_junction_support((replace(first, remapping_preset="sr"),), policy)
    assert mismatch.status is JunctionSupportStatus.DEFERRED
    assert "does not match" in mismatch.reasons[0]


def test_evidence_rejects_impossible_spanning_flank() -> None:
    with pytest.raises(InputValidationError, match="fit on both sides"):
        replace(evidence("case", 1), minimum_spanning_flank=51)


def test_reviewed_policy_requires_approved_review_artifact() -> None:
    with pytest.raises(InputValidationError, match="approved review artifact"):
        JunctionSupportPolicy("ont", "map-ont", 1, 0.0, 20, reviewed=True)
    with pytest.raises(InputValidationError, match="approved review artifact"):
        JunctionSupportPolicy(
            "ont",
            "map-ont",
            1,
            0.0,
            20,
            reviewed=True,
            review=replace(REVIEW, decision="rejected"),
        )


def test_reviewed_policy_matches_artifact_configuration_and_timestamp() -> None:
    with pytest.raises(InputValidationError, match="does not match"):
        JunctionSupportPolicy(
            "ont",
            "map-ont",
            4,
            0.3,
            20,
            reviewed=True,
            review=REVIEW,
            truth_metadata=METADATA,
            candidate_baseline_sha256="b" * 64,
        )
    with pytest.raises(InputValidationError, match="RFC 3339"):
        replace(REVIEW, reviewed_at="t")


def test_policy_review_loader_is_strict_and_typed(tmp_path: Path) -> None:
    path = tmp_path / "review.json"
    path.write_text(
        json.dumps(
            {
                "truth_dataset_sha256": "a" * 64,
                "candidate_baseline_sha256": "b" * 64,
                "reviewer": "reviewer",
                "reviewed_at": "2026-08-07T00:00:00Z",
                "decision": "approved",
                "technology": "ont",
                "remapping_preset": "map-ont",
                "minimum_spanning_reads": 3,
                "minimum_spanning_fraction": 0.3,
                "minimum_spanning_flank": 20,
                "minimum_mapping_quality": 0,
                "notes": "reviewed externally",
            }
        ),
        encoding="utf-8",
    )
    assert load_junction_policy_review(path) == replace(REVIEW, notes="reviewed externally")
    path.write_text("{}", encoding="utf-8")
    with pytest.raises(InputValidationError, match="missing"):
        load_junction_policy_review(path)
    invalid_type = {
        "truth_dataset_sha256": 1,
        "candidate_baseline_sha256": "b" * 64,
        "reviewer": "reviewer",
        "reviewed_at": "2026-08-07T00:00:00Z",
        "decision": "approved",
        "technology": "ont",
        "remapping_preset": "map-ont",
        "minimum_spanning_reads": 3,
        "minimum_spanning_fraction": 0.3,
        "minimum_spanning_flank": 20,
        "minimum_mapping_quality": 0,
    }
    path.write_text(json.dumps(invalid_type), encoding="utf-8")
    with pytest.raises(InputValidationError, match="must be a string"):
        load_junction_policy_review(path)
    path.write_text(
        '{"truth_dataset_sha256":"' + "a" * 64 + '","truth_dataset_sha256":"' + "a" * 64 + '"}',
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="duplicate JSON field"):
        load_junction_policy_review(path)


def test_truth_set_metadata_is_typed_and_unreviewed() -> None:
    metadata = load_junction_truth_set_metadata(
        Path(__file__).parents[1] / "benchmarks" / "pseudomonas_junction_truth_metadata.json",
        truth_path=DATASET / "expected" / "expected_junctions.tsv",
    )
    assert metadata.case_count == metadata.true_case_count + metadata.artificial_case_count
    assert metadata.false_support_baseline_established
    assert not metadata.reviewed


def test_policy_review_must_match_reviewed_truth_metadata() -> None:
    metadata = JunctionTruthSetMetadata(
        "dataset",
        "1",
        "ont",
        "map-ont",
        "reviewed truth",
        "a" * 64,
        19,
        15,
        4,
        True,
        reviewed=True,
    )
    review = replace(REVIEW, truth_dataset_sha256="a" * 64, candidate_baseline_sha256="b" * 64)
    validate_junction_policy_review(review, metadata, candidate_baseline_sha256="b" * 64)
    with pytest.raises(InputValidationError, match="truth digest"):
        validate_junction_policy_review(
            replace(review, truth_dataset_sha256="c" * 64),
            metadata,
            candidate_baseline_sha256="b" * 64,
        )
    with pytest.raises(InputValidationError, match="reviewed truth metadata"):
        validate_junction_policy_review(
            review, replace(metadata, reviewed=False), candidate_baseline_sha256="b" * 64
        )
