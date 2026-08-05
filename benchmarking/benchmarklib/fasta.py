import gzip

COMP = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


def reverse_complement(seq):
    return seq.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    records, name, chunks = [], None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError("sequence before FASTA header")
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def format_fasta(records, width=80):
    out = []
    for name, seq in records:
        out.append(f">{name}\n")
        out.extend(seq[i : i + width] + "\n" for i in range(0, len(seq), width))
    return "".join(out)
