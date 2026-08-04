"""Sequence transformations with explicit IUPAC handling."""

from contigger.exceptions import InputValidationError

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)
_IUPAC_DNA_MIXED_CASE = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an IUPAC DNA sequence."""
    invalid = sorted(set(sequence) - _IUPAC_DNA_MIXED_CASE)
    if invalid:
        raise InputValidationError("invalid IUPAC DNA symbol(s): " + ", ".join(invalid))
    return sequence.translate(_COMPLEMENT)[::-1]
