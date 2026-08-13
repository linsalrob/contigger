"""Validated resumable storage for complete classified relationship decisions."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from uuid import uuid4

from contigger.exceptions import InputValidationError
from contigger.models import (
    AlignmentHit,
    AlignmentType,
    CatalogueSequence,
    Orientation,
    PairRelationship,
    RejectedAlignment,
    Relationship,
    RelationshipType,
    RunConfig,
)

_FORMAT = 1


def relationship_artifact_path(output_prefix: Path) -> Path:
    """Return the persistent, non-biological relationship checkpoint path."""
    return output_prefix.parent / f".{output_prefix.name}.relationships.jsonl"


def relationship_artifact_identity(
    catalogue: Iterable[CatalogueSequence],
    config: RunConfig,
    *,
    minimap2_version: str | None = None,
) -> dict[str, object]:
    """Return the deterministic source/configuration identity for one checkpoint."""
    digest = hashlib.sha256()
    for sequence in sorted(catalogue, key=lambda item: item.identifier):
        digest.update(sequence.identifier.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sequence.sha256.encode("ascii"))
        digest.update(b"\0")
    return {
        "format": _FORMAT,
        "catalogue_sha256": digest.hexdigest(),
        "configuration": config.as_dict(),
        "minimap2_version": minimap2_version,
    }


def write_relationship_artifact(
    path: Path,
    identity: dict[str, object],
    relationships: Iterable[PairRelationship],
) -> None:
    """Atomically write sorted complete decisions for a later validated reuse."""
    with RelationshipArtifactWriter(path, identity) as output:
        for decision in sorted(
            relationships,
            key=lambda item: (item.relationship.target_id, item.relationship.query_id),
        ):
            output.write(decision)


class RelationshipArtifactWriter:
    """Atomically stream sorted pair decisions into a relationship checkpoint."""

    def __init__(self, path: Path, identity: dict[str, object]) -> None:
        """Prepare an artifact writer that records the supplied run identity."""
        self.path = path
        self.identity = identity
        self._temporary = path.parent / f".{path.name}.{uuid4().hex}.tmp"
        self._output = None
        self._last_key: tuple[str, str] | None = None
        self._decision_count = 0

    def __enter__(self) -> RelationshipArtifactWriter:
        """Open a temporary artifact and write its validation header."""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            self._output = self._temporary.open("w", encoding="utf-8", newline="")
            self._output.write(json.dumps({"identity": self.identity}, sort_keys=True) + "\n")
        except OSError as error:
            raise InputValidationError(
                f"cannot write relationship artifact {self.path}: {error}"
            ) from error
        return self

    def write(self, decision: PairRelationship) -> None:
        """Append one decision, requiring strictly increasing pair identifiers."""
        if self._output is None:
            raise InputValidationError("relationship artifact writer is not open")
        key = (decision.relationship.target_id, decision.relationship.query_id)
        if self._last_key is not None and key <= self._last_key:
            raise InputValidationError(
                "relationship artifact decisions must be written in unique sorted pair order"
            )
        self._output.write(json.dumps(_encode_pair(decision), sort_keys=True) + "\n")
        self._last_key = key
        self._decision_count += 1

    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> bool:
        """Publish only a complete artifact; remove a failed temporary artifact."""
        try:
            if self._output is not None:
                if exc_type is None:
                    self._output.write(
                        json.dumps(
                            {"complete": {"decision_count": self._decision_count}},
                            sort_keys=True,
                        )
                        + "\n"
                    )
                self._output.close()
            if exc_type is None:
                self._temporary.replace(self.path)
            elif self._temporary.exists():
                self._temporary.unlink()
        except OSError as error:
            raise InputValidationError(
                f"cannot finalize relationship artifact {self.path}: {error}"
            ) from error
        return False


def read_relationship_artifact(
    path: Path, identity: dict[str, object]
) -> tuple[PairRelationship, ...] | None:
    """Return a validated checkpoint, or ``None`` when it is absent or stale."""
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as input_file:
            header = json.loads(input_file.readline())
            if header != {"identity": identity}:
                return None
            payloads = tuple(json.loads(line) for line in input_file if line.strip())
    except (OSError, TypeError, ValueError, KeyError) as error:
        raise InputValidationError(f"cannot read relationship artifact {path}: {error}") from error
    try:
        if not payloads:
            raise ValueError("no completion record")
        footer = _mapping(payloads[-1])
        complete = footer.get("complete")
        if set(footer) != {"complete"} or not isinstance(complete, dict):
            raise ValueError("invalid completion record")
        expected_count = _integer(complete.get("decision_count"))
        if expected_count != len(payloads) - 1:
            raise ValueError("invalid decision count")
        decisions = tuple(_decode_pair(_mapping(payload)) for payload in payloads[:-1])
    except (TypeError, ValueError, KeyError) as error:
        raise InputValidationError(
            f"relationship artifact has an invalid completion record: {path}: {error}"
        ) from error
    keys = tuple((item.relationship.target_id, item.relationship.query_id) for item in decisions)
    if keys != tuple(sorted(keys)) or len(set(keys)) != len(keys):
        raise InputValidationError(
            f"relationship artifact is not deterministically ordered: {path}"
        )
    return tuple(
        sorted(
            decisions,
            key=lambda item: (item.relationship.query_id, item.relationship.target_id),
        )
    )


def _encode_pair(decision: PairRelationship) -> dict[str, object]:
    return {
        "relationship": _encode_relationship(decision.relationship),
        "representative_hit": _encode_hit(decision.representative_hit),
        "accepted_hits": [_encode_hit(item) for item in decision.accepted_hits],
        "rejected_alignments": [
            {"hit": _encode_hit(item.hit), "relationship": _encode_relationship(item.relationship)}
            for item in decision.rejected_alignments
        ],
        "ambiguity_reasons": list(decision.ambiguity_reasons),
    }


def _decode_pair(payload: dict[str, object]) -> PairRelationship:
    accepted = tuple(_decode_hit(item) for item in _items(payload, "accepted_hits"))
    rejected = tuple(
        RejectedAlignment(_decode_hit(item["hit"]), _decode_relationship(item["relationship"]))
        for item in _items(payload, "rejected_alignments")
    )
    representative = payload["representative_hit"]
    return PairRelationship(
        _decode_relationship(payload["relationship"]),
        None if representative is None else _decode_hit(representative),
        accepted,
        rejected,
        tuple(str(item) for item in _items(payload, "ambiguity_reasons")),
    )


def _encode_relationship(item: Relationship) -> dict[str, object]:
    return {
        "relationship_type": item.relationship_type.value,
        "query_id": item.query_id,
        "target_id": item.target_id,
        "orientation": item.orientation.value,
        "identity": item.identity,
        "aligned_length": item.aligned_length,
        "query_coverage": item.query_coverage,
        "target_coverage": item.target_coverage,
        "status": item.status,
        "reasons": list(item.reasons),
    }


def _decode_relationship(payload: object) -> Relationship:
    item = _mapping(payload)
    return Relationship(
        RelationshipType(_text(item["relationship_type"])),
        _text(item["query_id"]),
        _text(item["target_id"]),
        Orientation(_text(item["orientation"])),
        _number(item["identity"]),
        _integer(item["aligned_length"]),
        _number(item["query_coverage"]),
        _number(item["target_coverage"]),
        _text(item["status"]),
        tuple(_text(value) for value in _items(item, "reasons")),
    )


def _encode_hit(item: AlignmentHit | None) -> dict[str, object] | None:
    if item is None:
        return None
    return {
        "query_id": item.query_id,
        "target_id": item.target_id,
        "query_length": item.query_length,
        "target_length": item.target_length,
        "query_start": item.query_start,
        "query_end": item.query_end,
        "target_start": item.target_start,
        "target_end": item.target_end,
        "orientation": item.orientation.value,
        "matching_bases": item.matching_bases,
        "alignment_block_length": item.alignment_block_length,
        "mapping_quality": item.mapping_quality,
        "alignment_score": item.alignment_score,
        "supporting_seeds": item.supporting_seeds,
        "chaining_score": item.chaining_score,
        "secondary_chaining_score": item.secondary_chaining_score,
        "alignment_type": None if item.alignment_type is None else item.alignment_type.value,
    }


def _decode_hit(payload: object) -> AlignmentHit:
    item = _mapping(payload)
    alignment_type = item["alignment_type"]
    return AlignmentHit(
        _text(item["query_id"]),
        _text(item["target_id"]),
        _integer(item["query_length"]),
        _integer(item["target_length"]),
        _integer(item["query_start"]),
        _integer(item["query_end"]),
        _integer(item["target_start"]),
        _integer(item["target_end"]),
        Orientation(_text(item["orientation"])),
        _integer(item["matching_bases"]),
        _integer(item["alignment_block_length"]),
        _optional_int(item["mapping_quality"]),
        _optional_int(item["alignment_score"]),
        _optional_int(item["supporting_seeds"]),
        _optional_int(item["chaining_score"]),
        _optional_int(item["secondary_chaining_score"]),
        None if alignment_type is None else AlignmentType(_text(alignment_type)),
    )


def _mapping(value: object) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError("relationship artifact entry must be an object")
    return value


def _items(value: dict[str, object], key: str) -> list[object]:
    items = value[key]
    if not isinstance(items, list):
        raise ValueError(f"relationship artifact field {key!r} must be a list")
    return items


def _optional_int(value: object) -> int | None:
    return None if value is None else _integer(value)


def _integer(value: object) -> int:
    """Decode an integer JSON scalar without accepting arbitrary objects."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError("relationship artifact field must be an integer")
    return int(value)


def _number(value: object) -> float:
    """Decode a numeric JSON scalar without accepting arbitrary objects."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError("relationship artifact field must be numeric")
    return float(value)


def _text(value: object) -> str:
    """Decode a string JSON scalar without coercing structured values."""
    if not isinstance(value, str):
        raise ValueError("relationship artifact field must be text")
    return value
