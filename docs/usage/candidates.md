# `contigger candidates`

Emit positional-minimiser candidate evidence for inspection:

```bash
contigger candidates --manifest samples.tsv --output results/candidates.tsv
```

Candidate pairs are alignment requests, not relationships and not merges. Optional controls are `--kmer-size`, `--window-size`, `--min-shared-minimisers`, `--max-minimiser-frequency`, and `--terminal-band`. Use this command when investigating unexpectedly large or empty candidate sets.
