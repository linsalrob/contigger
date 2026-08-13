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
- [x] Record Slurm requested/used memory, elapsed time, CPU count, and exit reason.
- [x] Run a small fixture and one manageable Shark comparison to establish candidate,
      alignment, graph, and output counts.
- [x] Keep `--evidence none` for the initial scaling baseline so read validation does
      not add another source of cost.
- [x] Verify that exact/RC duplicates and containment remain unchanged as settings are
      tightened.

The first item is recorded in every successful `<output-prefix>.stats.json` file as
`input_contigs`, `input_bases`, `input_manifest`, and `input_by_sample`. The same file
also records peak process RSS, CPU seconds, and any Slurm allocation metadata available
to the process. Scheduler exit reasons for failed jobs remain an external Slurm
record; collect them with `sacct` using the job ID from the log filename.

For completed Setonix jobs, collect the authoritative scheduler record with:

```bash
scripts/collect_slurm_status.sh JOB_ID results/contigger.slurm.tsv
```

The TSV records Slurm state, exit code, elapsed time, allocated CPUs, maximum RSS,
requested memory, requested resources, and allocated resources. The collector must be
run after the job finishes because Slurm cannot report the final accounting state from
inside the running job.

### Milestone 1 completion evidence

The small fixture completed in sequence-only mode with 3 input contigs, 47 input bases,
zero candidate pairs, and 3 output contigs. A deterministic 4,000-contig subset of two
BundegiBeachWater assemblies was then run with minimap2 from the ATAVIDE environment:

| Setting | Candidates | Alignment observations | Graph components | Output contigs | Output bases |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline (`identity=98`, `min-overlap=1000`, `k=21`, `window=10`, `max-frequency=20`) | 47 | 47 | 4,000 | 4,000 | 2,108,711 |
| Exploratory tightened settings (`identity=99`, `min-overlap=2000`, `k=31`, `window=15`, `min-shared=8`, `max-frequency=20`) | 31 | 31 | 4,000 | 4,000 | 2,108,711 |

Both runs used `--evidence none`, completed successfully, constructed zero joins, and
produced byte-identical FASTA and provenance outputs. Catalogue and containment
dispositions were unchanged in both runs: 4,000 canonical sequences, 0 exact or
reverse-complement duplicate collapses, 1,927 reverse-oriented catalogue members, and
0 contained contigs removed.
The exploratory tightened run reduced candidate/alignment work on this subset, but it is
**not an accepted production setting**: the checked-in Pseudomonas decision-preservation
benchmark must pass before a threshold change can be adopted. The tightened evaluation
increases missed/false relationship decisions around the `identity_9799` regression, so
the baseline settings remain the only approved settings from this comparison.

The exact private-data manifest-construction command, selection rule, tool versions,
output checksums, and summary counts are recorded in
[`benchmarks/shark_baseline_results.json`](benchmarks/shark_baseline_results.json).
The manifest itself, source paths, and source/subset checksums are deliberately not
committed because the Shark input data are experimental.
This is a manageable baseline, not evidence that the original 500k-contig Shark
comparisons are now scalable; those remain Milestone 2 and 3 work.

## Milestone 2 — Prevent candidate explosion

- [x] Profile memory and runtime separately for catalogue loading and minimiser
      generation.
- [ ] Replace large Python minimiser/evidence object collections with compact integer
      IDs and bounded batches or shard files.
- [ ] Stream candidate generation instead of retaining every observation and pair in
      memory at once.
- [x] Add an early candidate-count/memory estimate and a clear configurable guardrail.
- [x] Benchmark `--max-minimiser-frequency`, k-mer size, and window size for recall
      versus candidate reduction; do not update biological baselines silently.

### Milestone 2 progress

The merge path now records deterministic minimiser-pressure counters in
`stats.json` under `candidate_generation`: all and retained minimiser observations,
discarded repetitive observations, unique minimisers, maximum per-pair evidence,
candidate pairs, and the frequency, retained-seed, pair-expansion, and candidate-filter
stage timings. Candidate generation uses a two-pass frequency/retained-seed workflow,
so it no longer keeps a full global observation tuple alongside the retained seed index.
Per-pair evidence is also accumulated as compact shared-value, position, orientation, and
count sets rather than retaining a tuple of full minimiser observations for every shared
seed. This preserves candidate geometry while removing the largest duplicated evidence
collection from ordinary runs.
The retained-seed index and pair-evidence map use compact integer sequence and minimiser
IDs internally; source identifiers and k-mer strings are retained only at the input/output
boundaries needed to preserve collision-safe candidate semantics.
`--max-candidate-pairs N` aborts before minimap2 alignment if the final candidate count
exceeds `N`, preventing accidental submission of an unbounded alignment job. This is an
initial operational guardrail, not yet a streaming candidate implementation; the
remaining unchecked tasks are still required to bound candidate-generation memory itself.
`--max-seed-pair-observations N` is an earlier guard: it rejects a run after the
frequency pass when the conservative per-minimiser seed-pair upper bound exceeds `N`,
before retained-seed indexing and pair expansion allocate their working structures.
Candidate-generation statistics now also report deterministic shard count plus temporary
seed and pair-evidence disk usage, so scale runs can quantify the memory/disk trade-off.
Completed `stats.json` files now also record current RSS before and after the catalogue
and candidate stages under `stage_resource_usage`, alongside their elapsed times. These
snapshots support profiling on Linux/Pawsey systems; they are not per-stage peak-memory
measurements and therefore do not replace Slurm `MaxRSS` collection.

Retained-seed shards are now externally sorted in fixed-size chunks and merged one
minimiser value at a time. This removes the previous whole-shard Python mapping of
retained observations while preserving deterministic candidate evidence. Temporary-sort
bytes are recorded separately. Final candidate tuples and per-pair evidence remain
materialized for the current public API, so this completes only the seed-group memory
bound; the remaining pair/graph streaming work belongs to Milestones 2–4.

### 10k sharded-candidate baseline (Setonix)

Job `46892769` completed on the `work` partition with 10,000 deterministic 100 bp
contigs (1,000,000 bases), 16 candidate shards, 8 requested CPUs, 16 GiB requested
memory, and no swaps. Candidate generation took 52.57 seconds (catalogue: 0.08 seconds)
and emitted 2,887 candidate pairs: a 17,317-fold reduction from all-vs-all pairs.
Maximum RSS was 723,672 KiB from `time -v` (Slurm `MaxRSS`: 857,504 KiB). Temporary seed
and pair shard files used 6,137,828 and 95,412,263 bytes respectively. Candidate filtering
was the dominant stage at 41.64 seconds; this is the main performance focus before a 100k
run. The complete deterministic result is tracked in
`benchmarks/scale-results/sharded-candidates-10k.json`.

### Pseudomonas candidate-parameter sweep (Setonix)

Job `46892891` completed with exit code `0:0`, no swaps, and 100,592 KiB maximum RSS
from `time -v` (Slurm `MaxRSS`: 155,436 KiB). It evaluated all 12 combinations of k-mer
sizes 21/31, minimiser windows 10/15, and maximum minimiser frequencies 20/50/100 against
both checked-in asm5 and asm20 PAF benchmarks. Every combination produced 61 candidate
pairs, complete candidate recall, two false relationship-stage merges, and respectively
four/two missed relationships. The small truth dataset therefore provides no evidence to
change defaults; the results are recorded in
`benchmarks/scale-results/pseudomonas-candidate-sweep.json`.

### Real-contig candidate scaling (Setonix)

The nested, identifier-hash-selected fixtures from a private 23,279,653-contig,
12,244,309,163-nucleotide assembly completed without swaps at 10k, 100k, and 1m
contigs. These runs measured only FASTA loading, catalogue construction, and candidate
generation (`k=21`, `window=10`, maximum minimiser frequency `20`, 64 shards, and a
100,000,000 potential seed-pair guard); they did not invoke minimap2 or construct
biological merges.

| Contigs | Bases | Candidates | Candidate time | Process peak RSS | Temporary seed/pair data |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10,000 | 5.23 Mb | 826 | 44.3 s | 257 MiB | 47.5 / 0.4 MB |
| 100,000 | 52.15 Mb | 21,884 | 460.8 s | 2.18 GiB | 480.8 / 13.1 MB |
| 1,000,000 | 524.48 Mb | 577,498 | 4,330.3 s | 20.74 GiB | 4.87 / 0.41 GB |

The 1m job completed with exit code `0:0` in 1h12m52s; Slurm reported 26,975,948 KiB
maximum RSS and approximately 5.03 GiB maximum disk write. Candidate time and retained
seed disk increased close to linearly, but candidate pairs and pair-evidence disk grew
faster than input count. The 1m fixture represents only 4.3% of the source contigs and
bases, so its runtime and memory must **not** be linearly extrapolated to the full
collection. The tracked aggregate result is
`benchmarks/scale-results/real-contig-candidate-scaling.json`; no private source path,
manifest, fixture, or raw scheduler log is committed.

### Private Shark manifest candidate-only preflights (Setonix)

The three smallest multi-assembly private manifests were profiled with the production
candidate implementation only: 64 candidate shards, a 100,000,000 seed-pair upper
bound, and a 1,000,000 final-candidate bound. No run invoked minimap2, built a graph, or
wrote a biological output. The aggregate, path-free results are tracked in
`benchmarks/shark_candidate_preflight_results.json`.

| Manifest | Input contigs | Candidates | Candidate time | Slurm MaxRSS | Outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| BundegiBeachWater | 502,796 | 567,392 | 42.9 min | 12.8 GiB | completed |
| Guitarfish | 775,558 | 546,286 | 50.6 min | 18.3 GiB | completed |
| TantabiddiWater | not recorded after guard rejection | 1,493,032 | 57.3 min | 19.6 GiB | rejected before alignment |

TantabiddiWater demonstrates why both guards must be explicit: its run stopped at the
final candidate-pair guard, before minimap2. BundegiBeachWater and Guitarfish passed
candidate generation but still have more than half a million candidate pairs each.
They are **not** approved for full merging until Milestones 3 and 4 provide bounded
alignment and graph processing.

## Milestone 3 — Make alignment batching scale

- [ ] Replace the current per-target index strategy for large runs with bounded
      multi-target index shards or another validated batched design.
- [ ] Preserve strict approved-pair filtering: alignments returned outside a planned
      batch must be rejected.
- [x] Reuse indexes only when sequence content, minimap2 preset, and software metadata
      match exactly.
- [x] Process alignment hits in bounded batches and write resumable relationship
      artifacts rather than retaining all hits in RAM.
- [ ] Report index builds, index reuse, alignment batches, temporary disk, and peak RSS.

### Milestone 3 progress

The indexed selective executor now splits each single-target approved-query group into
deterministic bounded batches (`--max-queries-per-alignment-batch`, default 1,000),
while reusing the same content-validated minimap2 index. Returned hits remain checked
against the exact approved query/target pairs in each batch. This bounds temporary query
FASTA size and individual minimap2 invocations without allowing multi-target or
all-vs-all expansion. Production merge now classifies each yielded batch immediately and
writes complete pair decisions to an atomic, digest/configuration-validated hidden
relationship artifact beside the output prefix. A compatible retry reuses that artifact
without launching minimap2, while stale artifacts are replaced and malformed or
truncated artifacts are rejected. Compatibility includes the exact minimap2 version as
well as catalogue content and normalized configuration. Raw alignment observations are
no longer accumulated globally. The current
graph and final output stages still load all classified relationships once, so incremental
component construction remains Milestone 4 work.

## Milestone 4 — Stream graph and sequence construction

- [ ] Build relationship components incrementally from relationship shards.
- [ ] Resolve and construct one bounded component at a time.
- [ ] Preserve deterministic output ordering and complete provenance across shards.
- [ ] Add safe resume points after candidate generation and alignment classification.
- [ ] Confirm that deferred, branched, cyclic, and conflict-containing components remain
      separate rather than being selected by score.

## Milestone 5 — Fix exact-match handling

- [ ] Reproduce the BundegiBeachWater failure on the private full manifest:
      `exact matches must be resolved by the sequence catalogue before graph construction`.
- [x] Ensure exact relationships are reconciled or filtered before graph construction.
- [x] Add a permanent synthetic regression test; rerun the exact/RC and Pseudomonas baselines.
- [ ] Confirm that the fix does not increase known false constructed joins.

## Milestone 6 — Validate at Shark scale

- [x] Run deterministic 10k and 100k-contig stress tests.
- [x] Extend the bounded real-contig candidate-stage stress test to 1m contigs and
      record RSS and temporary-disk pressure.
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
    --candidate-shards 64 \
    --max-seed-pair-observations 100000000 \
    --max-candidate-pairs 1000000 \
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
