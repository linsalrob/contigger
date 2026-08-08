# Performance and HPC

The production path reduces work in stages: exact/RC deduplication, positional-minimiser candidates, selective rather than all-vs-all alignment, and one-target indexed batches. Use `--index-dir` on fast scratch storage so validated target indexes can be reused within a run or across compatible reruns.

```bash
contigger merge --manifest samples.tsv --output-prefix results/contigger \
  --threads 32 --index-dir "$SCRATCH/contigger-indexes"
```

Example SLURM resources (not universal requirements):

```bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
contigger merge --manifest samples.tsv --output-prefix results/contigger \
  --threads "$SLURM_CPUS_PER_TASK"
```

Monitor `stats.json` for candidate pairs, alignment batches, index builds/reuse, stage timings, and output counts. The checked-in scale harness has been exercised at 10,000 contigs and a 100,000-contig smoke case; do not extrapolate those measurements to a different dataset or hardware without profiling. Unexpected candidate explosion is usually the first warning sign.
