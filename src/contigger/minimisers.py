"""Canonical positional-minimiser observations and conservative candidates."""

from __future__ import annotations

import hashlib
import time
from collections import Counter, deque
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TextIO

from contigger.exceptions import ConfigurationError, InputValidationError
from contigger.models import (
    CandidatePair,
    CatalogueSequence,
    MinimiserObservation,
    Orientation,
)
from contigger.utilities.sequences import reverse_complement

UNAMBIGUOUS_DNA = frozenset("ACGT")


@dataclass(frozen=True, slots=True)
class CandidateGenerationMetrics:
    """Candidate-generation pressure counters and stage timings."""

    input_sequences: int
    input_bases: int
    minimiser_observations: int
    retained_observations: int
    unique_minimisers: int
    repetitive_observations_discarded: int
    candidate_pairs: int
    maximum_pair_evidence: int
    frequency_pass_seconds: float
    retained_seed_pass_seconds: float
    pair_expansion_seconds: float
    candidate_filter_seconds: float

    def as_dict(self) -> dict[str, int | float]:
        """Return counters in the form used by the run statistics JSON."""
        return {
            "input_sequences": self.input_sequences,
            "input_bases": self.input_bases,
            "minimiser_observations": self.minimiser_observations,
            "retained_observations": self.retained_observations,
            "unique_minimisers": self.unique_minimisers,
            "repetitive_observations_discarded": self.repetitive_observations_discarded,
            "candidate_pairs": self.candidate_pairs,
            "maximum_pair_evidence": self.maximum_pair_evidence,
            "frequency_pass_seconds": self.frequency_pass_seconds,
            "retained_seed_pass_seconds": self.retained_seed_pass_seconds,
            "pair_expansion_seconds": self.pair_expansion_seconds,
            "candidate_filter_seconds": self.candidate_filter_seconds,
        }


def sequence_minimisers(
    sequence: CatalogueSequence, *, kmer_size: int, window_size: int
) -> tuple[MinimiserObservation, ...]:
    """Select every tied minimum in each window, retaining positions and strands."""
    _validate_parameters(kmer_size, window_size)
    if sequence.length < kmer_size:
        return ()
    kmer_count = sequence.length - kmer_size + 1
    effective_window = min(window_size, kmer_count)
    minima: deque[tuple[int, MinimiserObservation]] = deque()
    selected: set[MinimiserObservation] = set()
    for position in range(kmer_count):
        forward = sequence.sequence[position : position + kmer_size]
        if not set(forward) <= UNAMBIGUOUS_DNA:
            observation = None
        else:
            reverse = reverse_complement(forward)
            canonical = min(forward, reverse)
            orientation = Orientation.FORWARD if forward == canonical else Orientation.REVERSE
            value = int.from_bytes(hashlib.sha256(canonical.encode("ascii")).digest()[:8], "big")
            observation = MinimiserObservation(
                sequence_id=sequence.identifier,
                value=value,
                position=position,
                orientation=orientation,
                kmer=canonical,
            )

        if observation is not None:
            while minima and minima[-1][1].value > observation.value:
                minima.pop()
            minima.append((position, observation))
        window_start = position - effective_window + 1
        while minima and minima[0][0] < window_start:
            minima.popleft()
        if window_start < 0 or not minima:
            continue
        minimum_value = minima[0][1].value
        for _, minimum in minima:
            if minimum.value != minimum_value:
                break
            selected.add(minimum)
    return tuple(sorted(selected, key=_observation_sort_key))


def generate_candidates(
    sequences: Iterable[CatalogueSequence],
    *,
    kmer_size: int,
    window_size: int,
    min_shared_minimisers: int,
    max_minimiser_frequency: int,
    terminal_band: int,
) -> tuple[CandidatePair, ...]:
    """Generate candidates, retaining the historical tuple-only API."""
    candidates, _ = generate_candidates_with_metrics(
        sequences,
        kmer_size=kmer_size,
        window_size=window_size,
        min_shared_minimisers=min_shared_minimisers,
        max_minimiser_frequency=max_minimiser_frequency,
        terminal_band=terminal_band,
    )
    return candidates


def generate_candidates_with_metrics(
    sequences: Iterable[CatalogueSequence],
    *,
    kmer_size: int,
    window_size: int,
    min_shared_minimisers: int,
    max_minimiser_frequency: int,
    terminal_band: int,
) -> tuple[tuple[CandidatePair, ...], CandidateGenerationMetrics]:
    """Generate deterministic candidate pairs with explicit terminal seed geometry.

    Shared seeds are necessary but not sufficient: emitted pairs must also support a
    terminal overlap or end-to-end containment topology. Repetitive minimisers above
    the configured global observation frequency are discarded before pairing.
    """
    _validate_parameters(kmer_size, window_size)
    if min_shared_minimisers < 1 or max_minimiser_frequency < 1:
        raise ConfigurationError("minimiser frequency thresholds must be positive")
    if terminal_band < 0:
        raise ConfigurationError("terminal band cannot be negative")
    ordered = tuple(sorted(sequences, key=lambda item: item.identifier))
    if len({item.identifier for item in ordered}) != len(ordered):
        raise InputValidationError("candidate generation requires unique sequence identifiers")
    lengths = {item.identifier: item.length for item in ordered}
    frequency_started = time.monotonic()
    frequencies: Counter[tuple[int, str]] = Counter()
    minimiser_observations = 0
    for sequence in ordered:
        observations = sequence_minimisers(sequence, kmer_size=kmer_size, window_size=window_size)
        minimiser_observations += len(observations)
        frequencies.update((item.value, item.kmer) for item in observations)
    frequency_pass_seconds = time.monotonic() - frequency_started

    retained_started = time.monotonic()
    by_value: dict[tuple[int, str], list[MinimiserObservation]] = {}
    retained_observations = 0
    for sequence in ordered:
        for observation in sequence_minimisers(
            sequence, kmer_size=kmer_size, window_size=window_size
        ):
            key = (observation.value, observation.kmer)
            if frequencies[key] <= max_minimiser_frequency:
                by_value.setdefault(key, []).append(observation)
                retained_observations += 1
    retained_seed_pass_seconds = time.monotonic() - retained_started

    pair_expansion_started = time.monotonic()
    evidence: dict[tuple[str, str], list[tuple[MinimiserObservation, MinimiserObservation]]] = {}
    for value in sorted(by_value):
        items = sorted(by_value[value], key=_observation_sort_key)
        for left_index, left in enumerate(items):
            for right in items[left_index + 1 :]:
                if left.sequence_id == right.sequence_id:
                    continue
                query, target = sorted((left.sequence_id, right.sequence_id))
                first, second = (left, right) if left.sequence_id == query else (right, left)
                evidence.setdefault((query, target), []).append((first, second))
    pair_expansion_seconds = time.monotonic() - pair_expansion_started

    candidate_filter_started = time.monotonic()
    candidates: list[CandidatePair] = []
    for pair in sorted(evidence):
        matches = evidence[pair]
        shared_values = {(query.value, query.kmer) for query, _ in matches}
        if len(shared_values) < min_shared_minimisers:
            continue
        orientations = tuple(
            sorted(
                {
                    Orientation.FORWARD
                    if query.orientation is target.orientation
                    else Orientation.REVERSE
                    for query, target in matches
                },
                key=lambda item: item.value,
            )
        )
        topologies = _terminal_topologies(
            matches,
            lengths[pair[0]],
            lengths[pair[1]],
            kmer_size,
            terminal_band,
        )
        if not topologies:
            continue
        candidates.append(
            CandidatePair(
                query_id=pair[0],
                target_id=pair[1],
                shared_minimisers=len(shared_values),
                orientation=orientations[0] if len(orientations) == 1 else None,
                query_positions=tuple(sorted({query.position for query, _ in matches})),
                target_positions=tuple(sorted({target.position for _, target in matches})),
                supported_orientations=orientations,
                terminal_topologies=topologies,
                reasons=(
                    "canonical positional minimisers support terminal geometry",
                    *(() if len(orientations) == 1 else ("orientation evidence is ambiguous",)),
                ),
            )
        )
    ordered_candidates = tuple(candidates)
    candidate_filter_seconds = time.monotonic() - candidate_filter_started
    metrics = CandidateGenerationMetrics(
        input_sequences=len(ordered),
        input_bases=sum(item.length for item in ordered),
        minimiser_observations=minimiser_observations,
        retained_observations=retained_observations,
        unique_minimisers=len(frequencies),
        repetitive_observations_discarded=minimiser_observations - retained_observations,
        candidate_pairs=len(ordered_candidates),
        maximum_pair_evidence=max((len(matches) for matches in evidence.values()), default=0),
        frequency_pass_seconds=frequency_pass_seconds,
        retained_seed_pass_seconds=retained_seed_pass_seconds,
        pair_expansion_seconds=pair_expansion_seconds,
        candidate_filter_seconds=candidate_filter_seconds,
    )
    return ordered_candidates, metrics


def write_candidates_tsv(candidates: Iterable[CandidatePair], output: TextIO) -> None:
    """Write deterministic candidate evidence without implying a relationship."""
    columns = (
        "query",
        "target",
        "shared_minimisers",
        "orientation",
        "supported_orientations",
        "query_positions",
        "target_positions",
        "terminal_topologies",
        "reasons",
    )
    output.write("\t".join(columns) + "\n")
    for item in sorted(candidates, key=lambda candidate: (candidate.query_id, candidate.target_id)):
        output.write(
            "\t".join(
                (
                    item.query_id,
                    item.target_id,
                    str(item.shared_minimisers),
                    "." if item.orientation is None else item.orientation.value,
                    ",".join(value.value for value in item.supported_orientations),
                    ",".join(str(value) for value in item.query_positions),
                    ",".join(str(value) for value in item.target_positions),
                    ",".join(item.terminal_topologies),
                    "; ".join(item.reasons),
                )
            )
            + "\n"
        )


def _terminal_topologies(
    matches: list[tuple[MinimiserObservation, MinimiserObservation]],
    query_length: int,
    target_length: int,
    kmer_size: int,
    terminal_band: int,
) -> tuple[str, ...]:
    topologies: set[str] = set()

    def query_prefix(position: int) -> bool:
        return position < terminal_band

    def query_suffix(position: int) -> bool:
        return position + kmer_size > query_length - terminal_band

    def target_prefix(position: int) -> bool:
        return position < terminal_band

    def target_suffix(position: int) -> bool:
        return position + kmer_size > target_length - terminal_band

    by_orientation: dict[Orientation, list[tuple[MinimiserObservation, MinimiserObservation]]] = {
        Orientation.FORWARD: [],
        Orientation.REVERSE: [],
    }
    for query, target in matches:
        relative = (
            Orientation.FORWARD if query.orientation is target.orientation else Orientation.REVERSE
        )
        by_orientation[relative].append((query, target))

    for relative, oriented_matches in by_orientation.items():
        query_positions = {query.position for query, _ in oriented_matches}
        target_positions = {target.position for _, target in oriented_matches}
        if not query_positions:
            continue
        query_has_prefix = any(query_prefix(item) for item in query_positions)
        query_has_suffix = any(query_suffix(item) for item in query_positions)
        target_has_prefix = any(target_prefix(item) for item in target_positions)
        target_has_suffix = any(target_suffix(item) for item in target_positions)
        oriented_target_prefix = (
            target_has_prefix if relative is Orientation.FORWARD else target_has_suffix
        )
        oriented_target_suffix = (
            target_has_suffix if relative is Orientation.FORWARD else target_has_prefix
        )
        suffix = "" if relative is Orientation.FORWARD else "_REVERSE"
        if query_has_suffix and oriented_target_prefix:
            topologies.add(f"QUERY_SUFFIX_TO_TARGET_PREFIX{suffix}")
        if query_has_prefix and oriented_target_suffix:
            topologies.add(f"TARGET_SUFFIX_TO_QUERY_PREFIX{suffix}")
        if query_has_prefix and query_has_suffix:
            topologies.add(f"QUERY_CONTAINMENT_POSSIBLE{suffix}")
        if oriented_target_prefix and oriented_target_suffix:
            topologies.add(f"TARGET_CONTAINMENT_POSSIBLE{suffix}")
    return tuple(sorted(topologies))


def _validate_parameters(kmer_size: int, window_size: int) -> None:
    if kmer_size < 1 or window_size < 1:
        raise ConfigurationError("k-mer and minimiser window sizes must be positive")


def _observation_sort_key(item: MinimiserObservation) -> tuple[str, int, int, str, str]:
    return (item.sequence_id, item.position, item.value, item.orientation.value, item.kmer)
