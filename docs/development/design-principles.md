# Design principles

- A missed merge is preferable to a false merge.
- Identity alone never authorizes a join.
- Internal similarity is not terminal geometry.
- Orientation is explicit; coordinates are zero-based half-open.
- Pair ambiguity and graph ambiguity are preserved, never score-selected.
- Existing sample BAMs support source-contig observations, not a new junction by themselves.
- Contradictory sample-specific evidence is never pooled without a reviewed policy.
- No source contig disappears without a provenance disposition.
- Output ordering and IDs are deterministic.
- Any optimization must preserve candidate recall and final decisions on checked-in baselines.

The current safe construction layer emits only exact conflict-free terminal paths. Unsupported SNPs, indels, branches, repeats, cycles, and forbidden benchmark edges remain deferred.
