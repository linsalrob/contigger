# `contigger catalogue`

Create the exact/reverse-complement catalogue without aligning candidate pairs:

```bash
contigger catalogue --manifest samples.tsv \
  --output-fasta results/catalogue.fasta \
  --output-provenance results/catalogue.provenance.tsv
```

The catalogue writes one canonical sequence per exact strand-aware identity and records every source member. It does not construct overlaps or remove biological variation.
