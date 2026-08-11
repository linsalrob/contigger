"""Canonical positional-minimiser observations and conservative candidates."""

from __future__ import annotations

import hashlib
import time
from collections import Counter, deque
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
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
MAX_CANDIDATE_SHARDS = 64


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
    potential_seed_pair_observations: int
    frequency_pass_seconds: float
    retained_seed_pass_seconds: float
    pair_expansion_seconds: float
    candidate_filter_seconds: float
    candidate_shards: int
    temporary_seed_bytes: int
    temporary_pair_bytes: int

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
            "potential_seed_pair_observations": self.potential_seed_pair_observations,
            "frequency_pass_seconds": self.frequency_pass_seconds,
            "retained_seed_pass_seconds": self.retained_seed_pass_seconds,
            "pair_expansion_seconds": self.pair_expansion_seconds,
            "candidate_filter_seconds": self.candidate_filter_seconds,
            "candidate_shards": self.candidate_shards,
            "temporary_seed_bytes": self.temporary_seed_bytes,
            "temporary_pair_bytes": self.temporary_pair_bytes,
        }


@dataclass(slots=True)
class _PairEvidence:
    """Compact candidate evidence accumulated without retaining seed-pair tuples."""

    shared_values: set[int]
    query_positions: set[int]
    target_positions: set[int]
    positions_by_orientation: dict[Orientation, tuple[set[int], set[int]]]
    observation_count: int = 0

    @classmethod
    def empty(cls) -> _PairEvidence:
        """Create an empty accumulator for one ordered candidate pair."""
        return cls(set(), set(), set(), {})

    def add(
        self,
        minimiser_id: int,
        query_position: int,
        target_position: int,
        query_orientation: Orientation,
        target_orientation: Orientation,
    ) -> None:
        """Record one shared minimiser without retaining the source objects."""
        relative = (
            Orientation.FORWARD if query_orientation is target_orientation else Orientation.REVERSE
        )
        query_by_orientation, target_by_orientation = self.positions_by_orientation.setdefault(
            relative, (set(), set())
        )
        self.shared_values.add(minimiser_id)
        self.query_positions.add(query_position)
        self.target_positions.add(target_position)
        query_by_orientation.add(query_position)
        target_by_orientation.add(target_position)
        self.observation_count += 1


@dataclass(frozen=True, slots=True)
class _SeedObservation:
    """Retained minimiser evidence using a compact sequence index."""

    sequence_index: int
    position: int
    orientation: Orientation


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
    max_seed_pair_observations: int | None = None,
    candidate_shards: int = 16,
) -> tuple[CandidatePair, ...]:
    """Generate candidates, retaining the historical tuple-only API."""
    candidates, _ = generate_candidates_with_metrics(
        sequences,
        kmer_size=kmer_size,
        window_size=window_size,
        min_shared_minimisers=min_shared_minimisers,
        max_minimiser_frequency=max_minimiser_frequency,
        terminal_band=terminal_band,
        max_seed_pair_observations=max_seed_pair_observations,
        candidate_shards=candidate_shards,
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
    max_seed_pair_observations: int | None = None,
    candidate_shards: int = 16,
) -> tuple[tuple[CandidatePair, ...], CandidateGenerationMetrics]:
    """Generate deterministic candidate pairs with explicit terminal seed geometry.

    Shared seeds are necessary but not sufficient: emitted pairs must also support a
    terminal overlap or end-to-end containment topology. Repetitive minimisers above
    the configured global observation frequency are discarded before pairing.
    """
    _validate_parameters(kmer_size, window_size)
    if min_shared_minimisers < 1 or max_minimiser_frequency < 1:
        raise ConfigurationError("minimiser frequency thresholds must be positive")
    if max_seed_pair_observations is not None and max_seed_pair_observations < 1:
        raise ConfigurationError("maximum seed-pair observations must be positive when supplied")
    if candidate_shards < 1:
        raise ConfigurationError("candidate shard count must be positive")
    if candidate_shards > MAX_CANDIDATE_SHARDS:
        raise ConfigurationError(f"candidate shard count cannot exceed {MAX_CANDIDATE_SHARDS}")
    if terminal_band < 0:
        raise ConfigurationError("terminal band cannot be negative")
    ordered = tuple(sorted(sequences, key=lambda item: item.identifier))
    if len({item.identifier for item in ordered}) != len(ordered):
        raise InputValidationError("candidate generation requires unique sequence identifiers")
    sequence_ids = tuple(item.identifier for item in ordered)
    lengths = tuple(item.length for item in ordered)
    frequency_started = time.monotonic()
    frequencies: Counter[tuple[int, str]] = Counter()
    minimiser_observations = 0
    for sequence in ordered:
        observations = sequence_minimisers(sequence, kmer_size=kmer_size, window_size=window_size)
        minimiser_observations += len(observations)
        frequencies.update((item.value, item.kmer) for item in observations)
    frequency_pass_seconds = time.monotonic() - frequency_started
    potential_seed_pair_observations = sum(
        frequency * (frequency - 1) // 2
        for frequency in frequencies.values()
        if frequency <= max_minimiser_frequency
    )
    if (
        max_seed_pair_observations is not None
        and potential_seed_pair_observations > max_seed_pair_observations
    ):
        raise InputValidationError(
            "potential seed-pair observations "
            f"{potential_seed_pair_observations} exceed "
            f"--max-seed-pair-observations {max_seed_pair_observations}; "
            "tighten minimiser parameters or increase the limit"
        )

    retained_started = time.monotonic()
    retained_observations = 0
    temporary = TemporaryDirectory(prefix="contigger-candidates-")
    temporary_directory = temporary.name
    try:
        paths = [
            Path(temporary_directory) / f"seeds-{index:03d}.tsv"
            for index in range(candidate_shards)
        ]
        outputs = [path.open("w", encoding="ascii", newline="") for path in paths]
        try:
            for sequence_index, sequence in enumerate(ordered):
                for observation in sequence_minimisers(
                    sequence, kmer_size=kmer_size, window_size=window_size
                ):
                    key = (observation.value, observation.kmer)
                    if frequencies[key] <= max_minimiser_frequency:
                        shard = observation.value % candidate_shards
                        outputs[shard].write(
                            f"{observation.value}\t{observation.kmer}\t{sequence_index}\t"
                            f"{observation.position}\t{observation.orientation.value}\n"
                        )
                        retained_observations += 1
        finally:
            for output in outputs:
                output.close()
        retained_seed_pass_seconds = time.monotonic() - retained_started
        pair_expansion_started = time.monotonic()
        pair_paths = [
            Path(temporary_directory) / f"pairs-{index:03d}.tsv"
            for index in range(candidate_shards)
        ]
        pair_outputs = [path.open("w", encoding="ascii", newline="") for path in pair_paths]
        minimiser_id = 0
        try:
            for path in paths:
                by_value = _read_seed_shard(path)
                for value in sorted(by_value):
                    items = sorted(
                        by_value[value],
                        key=lambda item: (
                            item.sequence_index,
                            item.position,
                            item.orientation.value,
                        ),
                    )
                    for left_index, left in enumerate(items):
                        for right in items[left_index + 1 :]:
                            if left.sequence_index == right.sequence_index:
                                continue
                            query, target = sorted((left.sequence_index, right.sequence_index))
                            first, second = (
                                (left, right) if left.sequence_index == query else (right, left)
                            )
                            pair_shard = (query * len(sequence_ids) + target) % candidate_shards
                            pair_outputs[pair_shard].write(
                                f"{query}\t{target}\t{minimiser_id}\t{first.position}\t"
                                f"{second.position}\t{first.orientation.value}\t"
                                f"{second.orientation.value}\n"
                            )
                    minimiser_id += 1
        finally:
            for output in pair_outputs:
                output.close()
        pair_expansion_seconds = time.monotonic() - pair_expansion_started

        candidate_filter_started = time.monotonic()
        candidates: list[CandidatePair] = []
        maximum_pair_evidence = 0
        for path in pair_paths:
            evidence = _read_pair_shard(path)
            for pair in sorted(evidence):
                pair_evidence = evidence[pair]
                maximum_pair_evidence = max(maximum_pair_evidence, pair_evidence.observation_count)
                candidate = _candidate_from_evidence(
                    pair,
                    pair_evidence,
                    sequence_ids,
                    lengths,
                    kmer_size,
                    terminal_band,
                    min_shared_minimisers,
                )
                if candidate is not None:
                    candidates.append(candidate)
        temporary_seed_bytes = sum(path.stat().st_size for path in paths)
        temporary_pair_bytes = sum(path.stat().st_size for path in pair_paths)
    finally:
        temporary.cleanup()
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
        maximum_pair_evidence=maximum_pair_evidence,
        potential_seed_pair_observations=potential_seed_pair_observations,
        frequency_pass_seconds=frequency_pass_seconds,
        retained_seed_pass_seconds=retained_seed_pass_seconds,
        pair_expansion_seconds=pair_expansion_seconds,
        candidate_filter_seconds=candidate_filter_seconds,
        candidate_shards=candidate_shards,
        temporary_seed_bytes=temporary_seed_bytes,
        temporary_pair_bytes=temporary_pair_bytes,
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


def _read_seed_shard(path: Path) -> dict[tuple[int, str], list[_SeedObservation]]:
    """Read one deterministic retained-seed shard into its bounded grouping map."""
    by_value: dict[tuple[int, str], list[_SeedObservation]] = {}
    with path.open(encoding="ascii") as input_file:
        for line in input_file:
            value, kmer, sequence_index, position, orientation = line.rstrip("\n").split("\t")
            by_value.setdefault((int(value), kmer), []).append(
                _SeedObservation(int(sequence_index), int(position), Orientation(orientation))
            )
    return by_value


def _read_pair_shard(path: Path) -> dict[tuple[int, int], _PairEvidence]:
    """Reconstruct compact pair evidence from one bounded temporary shard."""
    evidence: dict[tuple[int, int], _PairEvidence] = {}
    with path.open(encoding="ascii") as input_file:
        for line in input_file:
            (
                query,
                target,
                minimiser,
                query_position,
                target_position,
                query_orientation,
                target_orientation,
            ) = line.rstrip("\n").split("\t")
            pair = (int(query), int(target))
            pair_evidence = evidence.get(pair)
            if pair_evidence is None:
                pair_evidence = _PairEvidence.empty()
                evidence[pair] = pair_evidence
            pair_evidence.add(
                int(minimiser),
                int(query_position),
                int(target_position),
                Orientation(query_orientation),
                Orientation(target_orientation),
            )
    return evidence


def _candidate_from_evidence(
    pair: tuple[int, int],
    evidence: _PairEvidence,
    sequence_ids: tuple[str, ...],
    lengths: tuple[int, ...],
    kmer_size: int,
    terminal_band: int,
    min_shared_minimisers: int,
) -> CandidatePair | None:
    """Build one candidate after its complete pair evidence has been reconstructed."""
    if len(evidence.shared_values) < min_shared_minimisers:
        return None
    orientations = tuple(sorted(evidence.positions_by_orientation, key=lambda item: item.value))
    topologies = _terminal_topologies_from_positions(
        evidence.positions_by_orientation,
        lengths[pair[0]],
        lengths[pair[1]],
        kmer_size,
        terminal_band,
    )
    if not topologies:
        return None
    return CandidatePair(
        query_id=sequence_ids[pair[0]],
        target_id=sequence_ids[pair[1]],
        shared_minimisers=len(evidence.shared_values),
        orientation=orientations[0] if len(orientations) == 1 else None,
        query_positions=tuple(sorted(evidence.query_positions)),
        target_positions=tuple(sorted(evidence.target_positions)),
        supported_orientations=orientations,
        terminal_topologies=topologies,
        reasons=(
            "canonical positional minimisers support terminal geometry",
            *(() if len(orientations) == 1 else ("orientation evidence is ambiguous",)),
        ),
    )


def _terminal_topologies(
    matches: list[tuple[MinimiserObservation, MinimiserObservation]],
    query_length: int,
    target_length: int,
    kmer_size: int,
    terminal_band: int,
) -> tuple[str, ...]:
    """Calculate topology labels from full observations for compatibility tests."""
    positions_by_orientation: dict[Orientation, tuple[set[int], set[int]]] = {}
    for query, target in matches:
        relative = (
            Orientation.FORWARD if query.orientation is target.orientation else Orientation.REVERSE
        )
        query_positions, target_positions = positions_by_orientation.setdefault(
            relative, (set(), set())
        )
        query_positions.add(query.position)
        target_positions.add(target.position)
    return _terminal_topologies_from_positions(
        positions_by_orientation,
        query_length,
        target_length,
        kmer_size,
        terminal_band,
    )


def _terminal_topologies_from_positions(
    positions_by_orientation: dict[Orientation, tuple[set[int], set[int]]],
    query_length: int,
    target_length: int,
    kmer_size: int,
    terminal_band: int,
) -> tuple[str, ...]:
    """Calculate terminal topology labels from compact per-orientation position sets."""
    topologies: set[str] = set()

    def query_prefix(position: int) -> bool:
        return position < terminal_band

    def query_suffix(position: int) -> bool:
        return position + kmer_size > query_length - terminal_band

    def target_prefix(position: int) -> bool:
        return position < terminal_band

    def target_suffix(position: int) -> bool:
        return position + kmer_size > target_length - terminal_band

    for relative, (query_positions, target_positions) in positions_by_orientation.items():
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
