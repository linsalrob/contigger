"""Small streaming-friendly FASTA parser used at validation boundaries."""

from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

from contigger.exceptions import FastaFormatError
from contigger.models import SequenceRecord

# DNA plus standard IUPAC ambiguity symbols. Gap characters are deliberately excluded.
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


def read_fasta(path: Path, sample: str = "") -> Iterator[SequenceRecord]:
    """Yield records from FASTA, uppercasing sequence and removing whitespace only.

    Whitespace anywhere on sequence lines is ignored. Every remaining character must
    be a standard IUPAC DNA symbol; characters are never silently discarded.
    """
    try:
        with path.open(encoding="utf-8") as handle:
            yield from _parse_fasta(handle, path, sample)
    except OSError as error:
        raise FastaFormatError(f"cannot read FASTA {path}: {error}") from error


def _parse_fasta(handle: TextIO, path: Path, sample: str) -> Iterator[SequenceRecord]:
    seen: set[str] = set()
    identifier: str | None = None
    description = ""
    chunks: list[str] = []
    ordinal = 0

    def make_record(line_number: int) -> SequenceRecord:
        nonlocal ordinal
        assert identifier is not None
        sequence = "".join(chunks).upper()
        if not sequence:
            raise FastaFormatError(f"{path}:{line_number}: empty sequence for {identifier!r}")
        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise FastaFormatError(
                f"{path}:{line_number}: invalid DNA symbol(s) for {identifier!r}: "
                + ", ".join(invalid)
            )
        record = SequenceRecord(
            identifier=f"{sample}:{identifier}" if sample else identifier,
            source_sample=sample,
            original_identifier=identifier,
            description=description,
            sequence=sequence,
            length=len(sequence),
            source_ordinal=ordinal,
        )
        ordinal += 1
        return record

    last_line = 0
    last_record_line = 0
    for line_number, raw_line in enumerate(handle, start=1):
        last_line = line_number
        line = raw_line.rstrip("\r\n")
        if line.startswith(">"):
            if identifier is not None:
                yield make_record(last_record_line)
            header = line[1:].strip()
            if not header:
                raise FastaFormatError(f"{path}:{line_number}: empty FASTA identifier")
            parts = header.split(maxsplit=1)
            identifier = parts[0]
            description = parts[1] if len(parts) == 2 else ""
            if identifier in seen:
                raise FastaFormatError(
                    f"{path}:{line_number}: duplicate FASTA identifier {identifier!r}"
                )
            seen.add(identifier)
            chunks = []
            last_record_line = line_number
        elif line.strip():
            if identifier is None:
                raise FastaFormatError(f"{path}:{line_number}: sequence data before first header")
            chunks.append("".join(line.split()))
            last_record_line = line_number
    if identifier is not None:
        yield make_record(last_record_line)
    elif last_line == 0:
        raise FastaFormatError(f"{path}: FASTA file is empty")


def validate_fasta(path: Path, sample: str = "") -> int:
    """Fully validate a FASTA file and return its record count."""
    return sum(1 for _ in read_fasta(path, sample))
