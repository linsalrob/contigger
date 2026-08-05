from .fasta import reverse_complement

def extract(seq, start, end, strand="+", circular=False):
    """Extract [start,end), allowing end > len(seq) only for circular references."""
    n=len(seq)
    if start < 0 or end < start: raise ValueError("invalid interval")
    if circular:
        if end-start > n: raise ValueError("interval longer than circle")
        part=(seq+seq)[start % n:start % n + end-start]
    else:
        if end > n: raise ValueError("interval outside reference")
        part=seq[start:end]
    if strand == "-": return reverse_complement(part)
    if strand != "+": raise ValueError("strand must be + or -")
    return part

def merge_suffix_prefix(left, right, overlap):
    if overlap < 0 or overlap > min(len(left),len(right)): raise ValueError("bad overlap")
    if left[-overlap:] != right[:overlap]: raise ValueError("overlap sequences differ")
    return left + right[overlap:]
