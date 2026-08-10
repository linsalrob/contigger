# Contigger scaling milestones

This checklist tracks the work needed to make large Contigger comparisons practical
while preserving the governing rule: **a missed merge is preferable to a false merge**.

The failed Shark runs showed three distinct classes of work: Python memory pressure,
alignment throughput, and one application-level exact-match failure. The missing
single-assembly manifests are intentionally excluded from this plan; those samples
do not need a cross-assembly comparison.

## Milestone 1 — Establish a reproducible baseline

- [x] Record the input contig count, total bases, and per-sample manifest for every
      comparison.
- [ ] Record Slurm requested/used memory, elapsed time, CPU count, and exit reason.
- [ ] Run a small fixture and one manageable Shark comparison to establish candidate,
      alignment, graph, and output counts.
- [ ] Keep `--evidence none` for the initial scaling baseline so read validation does
      not add another source of cost.
- [ ] Verify that exact/RC duplicates and containment remain unchanged as settings are
      tightened.

The first item is recorded in every successful `<output-prefix>.stats.json` file as
`input_contigs`, `input_bases`, `input_manifest`, and `input_by_sample`. The same file
also records peak process RSS, CPU seconds, and any Slurm allocation metadata available
to the process. Scheduler exit reasons for failed jobs remain an external Slurm
record; collect them with `sacct` using the job ID from the log filename.

## Milestone 2 — Prevent candidate explosion

- [ ] Profile memory and runtime separately for catalogue loading and minimiser
      generation.
- [ ] Replace large Python minimiser/evidence object collections with compact integer
      IDs and bounded batches or shard files.
- [ ] Stream candidate generation instead of retaining every observation and pair in
      memory at once.
- [ ] Add an early candidate-count/memory estimate and a clear configurable guardrail.
- [ ] Benchmark `--max-minimiser-frequency`, k-mer size, and window size for recall
      versus candidate reduction; do not update biological baselines silently.

## Milestone 3 — Make alignment batching scale

- [ ] Replace the current per-target index strategy for large runs with bounded
      multi-target index shards or another validated batched design.
- [ ] Preserve strict approved-pair filtering: alignments returned outside a planned
      batch must be rejected.
- [ ] Reuse indexes only when sequence content, minimap2 preset, and software metadata
      match exactly.
- [ ] Process alignment hits in bounded batches and write resumable relationship
      artifacts rather than retaining all hits in RAM.
- [ ] Report index builds, index reuse, alignment batches, temporary disk, and peak RSS.

## Milestone 4 — Stream graph and sequence construction

- [ ] Build relationship components incrementally from relationship shards.
- [ ] Resolve and construct one bounded component at a time.
- [ ] Preserve deterministic output ordering and complete provenance across shards.
- [ ] Add safe resume points after candidate generation and alignment classification.
- [ ] Confirm that deferred, branched, cyclic, and conflict-containing components remain
      separate rather than being selected by score.

## Milestone 5 — Fix exact-match handling

- [ ] Reproduce the BundegiBeachWater failure:
      `exact matches must be resolved by the sequence catalogue before graph construction`.
- [ ] Ensure exact relationships are reconciled or filtered before graph construction.
- [ ] Add a permanent regression test and rerun the exact/RC and Pseudomonas baselines.
- [ ] Confirm that the fix does not increase known false constructed joins.

## Milestone 6 — Validate at Shark scale

- [ ] Run deterministic 10k and 100k-contig stress tests.
- [ ] Run one representative multi-assembly Shark comparison under a documented Slurm
      allocation.
- [ ] Compare candidate counts, minimap2 process launches, index reuse, wall time, peak
      memory, temporary disk, merged paths, and deferred paths with the baseline.
- [ ] Increase scale only after the previous size completes without OOM or time limit.

## Recommended immediate operating settings

For an exploratory run, use a smaller manifest and conservative candidate controls,
then verify safe-merge recall:

```bash
contigger merge \
    --manifest samples.small.tsv \
    --output-prefix results/small \
    --identity 99 \
    --min-overlap 2000 \
    --kmer-size 31 \
    --window-size 15 \
    --min-shared-minimisers 8 \
    --max-minimiser-frequency 20 \
    --threads 16 \
    --evidence none \
    --index-dir "$MYSCRATCH/contigger-indexes/small"
```

These values are a diagnostic starting point, not universal biological defaults.
`--max-minimiser-frequency` is the most direct control over repetitive-seed pair
expansion. `--identity` and `--min-overlap` mainly affect later alignment and
relationship decisions.

## Definition of done

The scaling work is complete when a representative 100k-contig run finishes without
unbounded memory growth, uses bounded indexed alignment batches, produces deterministic
outputs with complete provenance, and retains zero false constructed joins on the
checked-in regressions. Only then should larger real metagenomic comparisons or broader
read-supported consensus policies be attempted.
