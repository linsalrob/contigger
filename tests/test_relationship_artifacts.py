"""Tests for atomic, validated relationship checkpoints."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from contigger.config import build_run_config
from contigger.exceptions import InputValidationError
from contigger.models import AlignmentHit, CatalogueSequence, Orientation, PairRelationship
from contigger.relationship_artifacts import (
    RelationshipArtifactWriter,
    read_relationship_artifact,
    relationship_artifact_identity,
    relationship_artifact_path,
    write_relationship_artifact,
)
from contigger.relationships import classify_pair


def _catalogue() -> tuple[CatalogueSequence, ...]:
    return (
        CatalogueSequence("a", "AACCGGTT", 8, "a" * 64, "a"),
        CatalogueSequence("b", "GGTTCCAA", 8, "b" * 64, "b"),
    )


def _decision() -> PairRelationship:
    config = build_run_config(min_overlap=4, min_containment=4, end_tolerance=0)
    hit = AlignmentHit("a", "b", 8, 8, 4, 8, 0, 4, Orientation.FORWARD, 4, 4)
    return classify_pair((hit,), config)


def test_relationship_artifact_round_trip_and_stale_identity(tmp_path: Path) -> None:
    """A checkpoint preserves complete pair decisions and rejects stale inputs."""
    config = build_run_config(output_prefix=tmp_path / "result")
    identity = relationship_artifact_identity(_catalogue(), config)
    path = relationship_artifact_path(config.output_prefix)
    decision = _decision()

    write_relationship_artifact(path, identity, (decision,))

    assert read_relationship_artifact(path, identity) == (decision,)
    changed = relationship_artifact_identity(
        _catalogue(), build_run_config(output_prefix=tmp_path / "other")
    )
    assert read_relationship_artifact(path, changed) is None


def test_relationship_artifact_writer_requires_sorted_unique_pairs(tmp_path: Path) -> None:
    """Streaming writers reject duplicate or out-of-order decisions before publish."""
    config = build_run_config(output_prefix=tmp_path / "result")
    path = relationship_artifact_path(config.output_prefix)
    with (
        pytest.raises(InputValidationError, match="unique sorted"),
        RelationshipArtifactWriter(
            path, relationship_artifact_identity(_catalogue(), config)
        ) as output,
    ):
        output.write(_decision())
        output.write(_decision())
    assert not path.exists()


def test_relationship_artifact_rejects_corrupt_content(tmp_path: Path) -> None:
    """A malformed matching checkpoint fails rather than being silently trusted."""
    config = build_run_config(output_prefix=tmp_path / "result")
    identity = relationship_artifact_identity(_catalogue(), config)
    path = relationship_artifact_path(config.output_prefix)
    path.write_text(json.dumps({"identity": identity}) + "\nnot-json\n", encoding="utf-8")

    with pytest.raises(InputValidationError, match="cannot read relationship artifact"):
        read_relationship_artifact(path, identity)


def test_relationship_artifact_rejects_a_validly_truncated_prefix(tmp_path: Path) -> None:
    """A missing completion record cannot silently drop a relationship decision."""
    config = build_run_config(output_prefix=tmp_path / "result")
    identity = relationship_artifact_identity(_catalogue(), config)
    path = relationship_artifact_path(config.output_prefix)
    write_relationship_artifact(path, identity, (_decision(),))
    path.write_text(
        "\n".join(path.read_text(encoding="utf-8").splitlines()[:-1]) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(InputValidationError, match="completion record"):
        read_relationship_artifact(path, identity)


def test_relationship_artifact_identity_includes_minimap2_version(tmp_path: Path) -> None:
    """A checkpoint produced by a different minimap2 release is stale."""
    config = build_run_config(output_prefix=tmp_path / "result")
    path = relationship_artifact_path(config.output_prefix)
    old_identity = relationship_artifact_identity(
        _catalogue(), config, minimap2_version="2.28-r1200"
    )
    write_relationship_artifact(path, old_identity, (_decision(),))

    new_identity = relationship_artifact_identity(
        _catalogue(), config, minimap2_version="2.29-r1300"
    )
    assert read_relationship_artifact(path, new_identity) is None
