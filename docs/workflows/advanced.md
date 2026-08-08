# Advanced workflow: evidence and scale

An advanced run may include multiple FASTA assemblies, indexed BAM/CRAM files, technology labels, assembly graphs, raw reads, and HPC storage. The merge command directly consumes the manifest FASTA and optional BAM/CRAM fields. It validates `technology` and `assembly_graph` paths for provenance, but does not use graph edges or raw FASTQ reads to authorize a merge.

Use explicit settings and fast scratch storage for indexes:

```bash
contigger merge --manifest samples.tsv --output-prefix results/contigger \
  --identity 98 --min-overlap 1000 --threads 32 \
  --minimap2-preset asm20 --index-dir "$SCRATCH/contigger-indexes" \
  --evidence alignments --emit-gfa
```

`asm5`, `asm10`, and `asm20` are passed to minimap2; the default is `asm20`. Index metadata records target content, preset, and minimap2 version. Stale or mismatched indexes are rejected rather than reused. `stats.json` records index builds, reuse, alignment batches, tool versions, and stage timings.

Keep raw reads and source BAMs separately. The current `benchmark-junction-remapping` command is a checked-in benchmark workflow, not a general raw-read merge input. No current reviewed policy turns its observations into merge authorization. For large jobs, monitor candidate counts and peak memory, keep originals, and inspect deferred components before changing thresholds.
