import random

BASES="ACGT"
def substitution_positions(seq, count, seed):
    eligible=[i for i,b in enumerate(seq.upper()) if b in BASES]
    if count > len(eligible): raise ValueError("too many substitutions")
    return sorted(random.Random(seed).sample(eligible,count))

def substitute(seq, count, seed):
    chars=list(seq.upper()); rows=[]
    rng=random.Random(seed)
    pos=substitution_positions(seq,count,seed)
    for p in pos:
        old=chars[p]; choices=[b for b in BASES if b != old]
        new=choices[rng.randrange(3)]; chars[p]=new; rows.append((p,old,new))
    return "".join(chars), rows

def insert(seq, position, bases): return seq[:position]+bases+seq[position:]
def delete(seq, position, length): return seq[:position]+seq[position+length:]
