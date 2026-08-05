import hashlib


def stable_hash(name, seed=47291):
    return int.from_bytes(hashlib.sha256(f"{seed}\0{name}".encode()).digest()[:8], "big")


def partition(name, partitions=3, seed=47291):
    return stable_hash(name, seed) % partitions


def read_fastq(handle):
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        plus = handle.readline()
        q = handle.readline()
        if not (s and plus and q) or not h.startswith("@") or not plus.startswith("+"):
            raise ValueError("malformed FASTQ")
        s = s.rstrip("\r\n")
        q = q.rstrip("\r\n")
        if len(s) != len(q):
            raise ValueError("FASTQ sequence/quality length mismatch")
        yield h[1:].split()[0], s, q
