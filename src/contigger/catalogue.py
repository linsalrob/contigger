"""Stable exact and reverse-complement sequence catalogue construction."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from pathlib import Path
from typing import TextIO

from contigger.exceptions import InputValidationError
from contigger.fasta import read_fasta
from contigger.models import (
    CatalogueMember,
    CatalogueSequence,
    Orientation,
    SampleInput,
    SequenceCatalogue,
    SequenceRecord,
)
from contigger.provenance import ProvenanceRecord
from contigger.utilities.sequences import reverse_complement


def load_source_sequences(samples: Iterable[SampleInput]) -> tuple[SequenceRecord, ...]:
    """Load all sample FASTAs with globally unique sample-qualified identifiers."""
    records: list[SequenceRecord] = []
    seen: set[str] = set()
    for sample in sorted(samples, key=lambda item: item.sample):
        for record in read_fasta(sample.contigs, sample.sample):
            if record.identifier in seen:
                raise InputValidationError(
                    f"duplicate sample-qualified source identifier: {record.identifier}"
                )
            seen.add(record.identifier)
            records.append(record)
    return tuple(sorted(records, key=_source_sort_key))


def build_catalogue(records: Iterable[SequenceRecord]) -> SequenceCatalogue:
    """Collapse only byte-exact forward or reverse-complement sequence duplicates."""
    ordered = tuple(sorted(records, key=_source_sort_key))
    source_ids = [record.identifier for record in ordered]
    if len(source_ids) != len(set(source_ids)):
        duplicate = next(item for item in source_ids if source_ids.count(item) > 1)
        raise InputValidationError(f"duplicate source identifier in catalogue input: {duplicate}")

    groups: dict[str, list[tuple[SequenceRecord, Orientation]]] = {}
    canonical_sequences: dict[str, str] = {}
    for record in ordered:
        reverse = reverse_complement(record.sequence)
        canonical = min(record.sequence, reverse)
        orientation = Orientation.FORWARD if record.sequence == canonical else Orientation.REVERSE
        digest = hashlib.sha256(canonical.encode("ascii")).hexdigest()
        previous = canonical_sequences.setdefault(digest, canonical)
        if previous != canonical:
            raise InputValidationError(f"SHA-256 collision while cataloguing {record.identifier}")
        groups.setdefault(digest, []).append((record, orientation))

    sequences: list[CatalogueSequence] = []
    members: list[CatalogueMember] = []
    for digest in sorted(groups):
        catalogue_id = f"ctg_{digest}"
        group = sorted(groups[digest], key=lambda item: _source_sort_key(item[0]))
        representative = min(
            group,
            key=lambda item: (
                item[1] is Orientation.REVERSE,
                *_source_sort_key(item[0]),
            ),
        )[0]
        sequence = canonical_sequences[digest]
        sequences.append(
            CatalogueSequence(
                identifier=catalogue_id,
                sequence=sequence,
                length=len(sequence),
                sha256=digest,
                representative_source_id=representative.identifier,
            )
        )
        for record, orientation in group:
            members.append(
                CatalogueMember(
                    catalogue_id=catalogue_id,
                    source_id=record.identifier,
                    source_sample=record.source_sample,
                    original_identifier=record.original_identifier,
                    orientation=orientation,
                    representative=record.identifier == representative.identifier,
                )
            )
    return SequenceCatalogue(
        sequences=tuple(sorted(sequences, key=lambda item: item.identifier)),
        members=tuple(sorted(members, key=lambda item: (item.catalogue_id, item.source_id))),
    )


def catalogue_provenance(catalogue: SequenceCatalogue) -> tuple[ProvenanceRecord, ...]:
    """Return one full-span provenance record for every source catalogue member."""
    lengths = {sequence.identifier: sequence.length for sequence in catalogue.sequences}
    group_sizes: dict[str, int] = {}
    for member in catalogue.members:
        group_sizes[member.catalogue_id] = group_sizes.get(member.catalogue_id, 0) + 1
    return tuple(
        ProvenanceRecord(
            output_sequence=member.catalogue_id,
            source_sample=member.source_sample,
            source_contig=member.original_identifier,
            relationship="EXACT_MATCH" if group_sizes[member.catalogue_id] > 1 else "UNIQUE",
            orientation=member.orientation,
            source_start=0,
            source_end=lengths[member.catalogue_id],
            output_start=0,
            output_end=lengths[member.catalogue_id],
            identity=1.0,
            disposition="representative" if member.representative else "exact_duplicate",
            decision_reason="canonical exact sequence catalogue",
        )
        for member in catalogue.members
    )


def write_catalogue_fasta(catalogue: SequenceCatalogue, output: TextIO) -> None:
    """Write canonical catalogue sequences in deterministic FASTA order."""
    for sequence in catalogue.sequences:
        output.write(f">{sequence.identifier}\n")
        for start in range(0, sequence.length, 80):
            output.write(sequence.sequence[start : start + 80] + "\n")


def write_catalogue_fasta_path(catalogue: SequenceCatalogue, path: Path) -> None:
    """Write a catalogue FASTA path with a domain-specific failure."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="ascii", newline="") as output:
            write_catalogue_fasta(catalogue, output)
    except OSError as error:
        raise InputValidationError(f"cannot write catalogue FASTA {path}: {error}") from error


def _source_sort_key(record: SequenceRecord) -> tuple[str, str, int, str]:
    return (
        record.source_sample,
        record.original_identifier,
        -1 if record.source_ordinal is None else record.source_ordinal,
        record.identifier,
    )
