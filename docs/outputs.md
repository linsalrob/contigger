# Outputs

For prefix `results/contigger`, the merge writes:

| File | Purpose |
| --- | --- |
| `.fasta` | Exact/RC representatives, eligible containment representatives, safe constructed paths, and retained deferred contigs. |
| `.provenance.tsv` | One or more rows tracing every source contig to output coordinates and disposition. |
| `.relationships.tsv` | Pair classifications, orientation, identity, length, status, and reasons. |
| `.ambiguous.tsv` | Deferred components, paths, conflicts, cycles, and reasons. |
| `.gfa` | GFA header, emitted segments, and eligible links when `--emit-gfa` is used; otherwise an empty file. |
| `.stats.json` | Counts, configuration, tool versions/commands, index metrics, and stage timings. |
| `.join_support.tsv` | Alignment-evidence join diagnostics; currently imperfect joins are reported as deferred. |
| `.variants.tsv` | Variant-resolution schema; currently empty because no unreviewed consensus is invented. |

The provenance table includes `output_id`, source sample/contig, orientation, zero-based half-open source and output coordinates, disposition, relationship, identity, reason, and evidence mode. No source contig should disappear without a provenance disposition.

An ambiguous row is useful information: it means Contigger deliberately retained alternatives rather than selecting a branch by score.
