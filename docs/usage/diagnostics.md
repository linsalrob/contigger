# Diagnostics and benchmarks

The following commands are primarily for reproducibility and developers:

```bash
contigger classify-paf --paf observations.paf.gz --output relationships.tsv
contigger benchmark --dataset test_data --paf test_data/alignments/all_vs_all.asm20.paf.gz
contigger benchmark-pipeline --dataset test_data --paf test_data/alignments/all_vs_all.asm20.paf.gz
contigger benchmark-junction-remapping --dataset test_data --preset map-ont
```

`classify-paf` classifies complete ordered pairs; it does not merge sequences. Benchmark commands compare deterministic checked-in observations and should not be treated as a substitute for inspecting a real merge's provenance.
