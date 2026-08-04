"""Deterministic small synthetic contigs for relationship evaluation only."""

from __future__ import annotations

import random
from dataclasses import dataclass

from contigger.models import RelationshipType


@dataclass(frozen=True, slots=True)
class SyntheticCase:
    """One independently alignable synthetic query-target truth case."""

    name: str
    query: str
    target: str
    truth: RelationshipType


def synthetic_cases(seed: int = 1729) -> tuple[SyntheticCase, ...]:
    """Return reproducible fixtures spanning conservative boundary conditions."""
    rng = random.Random(seed)
    core = _dna(rng, 400)
    left = _dna(rng, 250)
    right = _dna(rng, 250)
    overlap = _dna(rng, 300)
    changed = list(overlap)
    for position in range(0, len(changed), 50):
        changed[position] = _different_base(changed[position])
    indel_overlap = overlap[:150] + "A" + overlap[150:]
    repeat = _dna(rng, 140)
    repeat_query = repeat + _dna(rng, 180) + repeat
    repeat_target = repeat + _dna(rng, 220) + repeat
    internal = _dna(rng, 170)
    return (
        SyntheticCase("exact_duplicate", core, core, RelationshipType.EXACT_MATCH),
        SyntheticCase(
            "reverse_complement_duplicate",
            _reverse_complement(core),
            core,
            RelationshipType.EXACT_MATCH,
        ),
        SyntheticCase(
            "query_contained", core, left + core + right, RelationshipType.QUERY_CONTAINED_IN_TARGET
        ),
        SyntheticCase(
            "target_contained",
            left + core + right,
            core,
            RelationshipType.TARGET_CONTAINED_IN_QUERY,
        ),
        SyntheticCase(
            "forward_suffix_prefix",
            left + overlap,
            overlap + right,
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        ),
        SyntheticCase(
            "reverse_suffix_prefix",
            left + _reverse_complement(overlap),
            right + overlap,
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        ),
        SyntheticCase(
            "identity_98_percent",
            left + "".join(changed),
            overlap + right,
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        ),
        SyntheticCase(
            "small_indel",
            left + indel_overlap,
            overlap + right,
            RelationshipType.QUERY_SUFFIX_TO_TARGET_PREFIX,
        ),
        SyntheticCase(
            "internal_conserved",
            left + internal + right,
            _dna(rng, 200) + internal + _dna(rng, 200),
            RelationshipType.NO_RELATIONSHIP,
        ),
        SyntheticCase(
            "repeated_hits", repeat_query, repeat_target, RelationshipType.AMBIGUOUS_OVERLAP
        ),
        SyntheticCase(
            "incompatible_terminal_placements",
            repeat + _dna(rng, 160) + repeat,
            repeat + _dna(rng, 160) + repeat,
            RelationshipType.AMBIGUOUS_OVERLAP,
        ),
        SyntheticCase("low_complexity", "A" * 400, "C" * 400, RelationshipType.NO_RELATIONSHIP),
        SyntheticCase(
            "unrelated", _dna(rng, 400), _dna(rng, 400), RelationshipType.NO_RELATIONSHIP
        ),
    )


def mutate_substitutions(sequence: str, count: int) -> str:
    """Mutate the first ``count`` positions deterministically for boundary fixtures."""
    bases = list(sequence)
    for position in range(count):
        bases[position] = _different_base(bases[position])
    return "".join(bases)


def _dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _different_base(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
