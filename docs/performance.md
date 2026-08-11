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

For an unfamiliar collection, use a candidate guardrail based on a small representative
run. `--max-seed-pair-observations N` stops after frequency counting but before retained
seed indexing if a conservative seed-pair upper bound exceeds `N`; it protects the
candidate-generation working set. `--max-candidate-pairs N` stops before minimap2
alignment if more than `N` pairs survive minimiser filtering. Completed runs report `candidate_generation` counters in
`stats.json`, including total and retained minimiser observations, repetitive seeds
discarded, maximum evidence accumulated for one pair, and timings for the minimiser
frequency, retained-seed, pair-expansion, and candidate-filter stages. The current
two-pass implementation avoids retaining a second full global observation collection,
and stores compact per-pair seed summaries instead of full observation-pair tuples,
but it does not yet make candidate generation fully streaming; the guard remains an
alignment-cost safety valve rather than a substitute for the remaining scaling work.

For profiling, `stats.json` also records elapsed stage timings plus current RSS before
and after the catalogue and candidate stages in `stage_resource_usage`. These are useful
process snapshots on Linux, not replacement per-stage peak measurements; retain the
Slurm accounting record for final `MaxRSS`.

After a Slurm job completes, collect scheduler-side memory, elapsed-time, CPU, and
exit-status data separately from Contigger's process stats:

```bash
scripts/collect_slurm_status.sh JOB_ID results/contigger.slurm.tsv
```

This uses `sacct` and should be run after completion so `MaxRSS` and the final state are
available. Failed jobs may not produce a final Contigger `stats.json`, but this TSV
still records their scheduler exit code and state.

## Representative real-data scale fixtures

Do not test a large collection by copying it into the repository.  The optional
`benchmarking/real_fasta_scale.py` helper makes deterministic nested FASTA subsets from
one source FASTA, retaining the original contig sequences and identifiers.  Its
identifier-hash selection is independent of input order and approximately preserves the
source contig-length distribution.  For example, create 10,000-, 100,000-, and
1,000,000-contig fixtures on scratch storage:

```bash
python benchmarking/real_fasta_scale.py \
  --source /path/to/final.contigs.fa \
  --output-directory "$MYSCRATCH/contigger-scale-fixtures" \
  --records 10000 100000 1000000 \
  --seed contigger-scale-v1
```

The helper writes a manifest for each fixture plus `real-contigs-fixtures.json`, which
records the measured source record/base totals, seed, selection method, and fixture
sizes.  Keep these derived files and their Slurm accounting records with the private
dataset; do not commit them.  Profile the actual catalogue and candidate implementation
before running minimap2:

```bash
python benchmarking/profile_fasta_candidates.py \
  --fasta "$MYSCRATCH/contigger-scale-fixtures/real-contigs-100000.fasta" \
  --output results/real-candidates-100000.json \
  --max-minimiser-frequency 20 \
  --candidate-shards 64 \
  --max-seed-pair-observations 100000000
```

Run each size under a bounded Slurm allocation.  Increase the record count only after
checking candidate count, temporary shard bytes, peak RSS, and elapsed time at the
previous size.  A hash sample measures typical sequence-length and repeat content; it
does not guarantee a worst-case repeat stress test.

One private metagenomic assembly was profiled with this method at 10k, 100k, and 1m
contigs. The candidate stage took 44 seconds, 7.7 minutes, and 72.2 minutes,
respectively; process peak RSS was 257 MiB, 2.18 GiB, and 20.74 GiB. Candidate pairs
grew from 826 to 21,884 to 577,498, so pair evidence—not just input bases—must remain a
hard operational guardrail. Across these three nested fixtures, candidate time increased
10.4x and then 9.4x for each 10x increase in contig count; peak process RSS increased
8.7x and then 9.5x. The current retained-seed implementation therefore shows
approximately linear time and memory behaviour on this real dataset through 1m contigs.
This is an observation, not a general complexity guarantee: repetitive sequence content
can still make candidate-pair growth super-linear. The complete aggregate settings and
measurements are in `benchmarks/scale-results/real-contig-candidate-scaling.json`.
This was a pre-alignment benchmark: it does not include minimap2, graph, or output costs
and must not be extrapolated linearly to a whole collection.

For a Pawsey Setonix CPU example, see [`scripts/run_contigger_setonix.slurm`](https://github.com/linsalrob/contigger/blob/main/scripts/run_contigger_setonix.slurm). Replace its account placeholder before submission. It follows Pawsey's guidance to specify nodes, tasks, CPUs, memory, and wall time; see the [Setonix job documentation](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix) for current partition and allocation policies.
