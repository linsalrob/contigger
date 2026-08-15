# `contigger merge`

Merge is the main biological command. It validates the manifest, exact/RC-catalogues sequences, generates positional-minimiser candidates, aligns only planned pairs, classifies relationships, applies conservative graph decisions, and writes outputs.

```bash
contigger merge --manifest samples.tsv --output-prefix results/contigger \
  --identity 98 --evidence none --threads 16
```

Important options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--identity` | 98 | Minimum percent identity for relationship classification; not merge permission. |
| `--min-overlap` | 1000 | Minimum terminal overlap length. |
| `--min-containment` | 500 | Minimum containment span. |
| `--containment-coverage` | 98 | Required percent coverage of the contained sequence. |
| `--end-tolerance` | 50 | Allowed terminal coordinate slack. |
| `--threads` | 1 | minimap2 threads. |
| `--minimap2-preset` | `asm20` | `asm5`, `asm10`, or `asm20`. |
| `--index-dir` | output-local cache | Directory for validated target indexes. |
| `--max-candidate-pairs` | unlimited | Abort before alignment if this many candidate pairs is exceeded. |
| `--candidate-shards` | 16 | Disk shards used for candidate evidence; use 64 for the bounded large-scale profile. |
| `--max-queries-per-alignment-batch` | 1000 | Maximum approved queries passed to one indexed minimap2 invocation. |
| `--evidence` | `none` | `none`, `alignments`, or legacy `reads` (the latter is rejected). |
| `--emit-gfa` | off | Populate the GFA output. |
| `--dry-run` | off | Validate and print the plan without biological outputs. |

`--evidence none` permits exact/RC duplicates, eligible containment, and exact conflict-free terminal paths. `--evidence alignments` requires a validated BAM/CRAM for every sample and records imperfect-junction diagnostics; no currently reviewed policy authorizes a consensus base, so SNP/indel joins remain deferred. Conflict-policy values are accepted for compatibility but do not override conservative graph decisions.

During a normal run, Contigger prints five stage summaries—catalogue, candidates,
alignment/classification, graph planning, and construction/output preparation. Indexed
alignment also reports completed batches at the first batch, five-percent increments,
30-second intervals, and completion. These lines are intentionally throttled so large
runs remain readable while still showing whether work is advancing.

### Relationship checkpoints

During alignment, Contigger writes complete classified pair decisions to the atomic,
hidden checkpoint beside the output prefix: `.contigger.relationships.jsonl` for an
`--output-prefix results/contigger`. A later run with the same canonical input content
and normalized configuration reuses this checkpoint instead of rerunning minimap2.
The checkpoint is validated against the catalogue digests, normalized configuration, and
minimap2 version, and carries a completion record so a truncated prefix cannot be reused.
An old checkpoint is automatically replaced, while a malformed checkpoint stops the run
rather than being trusted. Keep it with the index directory when a large run may need restarting.
Delete it deliberately only when you want to force relationship classification again.
