# Parameters

Defaults favor fewer false joins. Increase thresholds to become more conservative; decrease them only with benchmark evidence.

| Group | Parameters | Guidance |
| --- | --- | --- |
| Relationship | `--identity`, `--min-overlap`, `--min-containment`, `--containment-coverage`, `--end-tolerance` | Higher identity/length/coverage or lower tolerance generally reduces false positives and misses more joins. |
| Candidates | `--kmer-size`, `--window-size`, `--min-shared-minimisers`, `--max-minimiser-frequency`, `--max-candidate-pairs` | Candidate settings control recall and runtime, not merge authorization. Leave defaults alone initially. |
| Performance | `--threads`, `--minimap2-preset`, `--index-dir` | Threads affect alignment runtime; presets affect sensitivity. Index metadata prevents unsafe reuse. |
| Evidence | `--evidence`, `--conflict-policy` | Alignment mode validates sample evidence but currently defers unreviewed imperfect consensus. |

Defaults are identity 98, overlap 1000, containment 500, containment coverage 98, end tolerance 50, k-mer 21, window 10, five shared minimisers, maximum minimiser frequency 100, one thread, and preset `asm20`.

`--max-candidate-pairs` has no default limit. Set it to a positive number on large or
unfamiliar collections to stop the run before minimap2 alignment when the candidate
set is unexpectedly large. The completed-run `stats.json` records minimiser observation
and candidate-pressure counters to help choose a safe value for later runs.
