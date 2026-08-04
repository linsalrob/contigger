"""Sequence transformations with explicit IUPAC handling."""

from contigger.exceptions import InputValidationError

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an IUPAC DNA sequence."""
    valid = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
    invalid = sorted(set(sequence) - valid)
    if invalid:
        raise InputValidationError("invalid IUPAC DNA symbol(s): " + ", ".join(invalid))
    return sequence.translate(_COMPLEMENT)[::-1]
