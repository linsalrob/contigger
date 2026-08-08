# Architecture

The production pipeline is:

```text
manifest → source validation → catalogue → exact/RC deduplication
→ positional minimiser candidates → indexed selective minimap2 alignment
→ complete-pair relationship classification → ambiguity-preserving graph
→ decision policy → provenance-complete path planning → sequence construction
→ FASTA, provenance, relationships, ambiguity, GFA, stats, evidence diagnostics
```

`catalogue.py` provides canonical sequences and source membership. `minimisers.py` emits positional candidate evidence. `alignment_planning.py` groups approved queries by exactly one target and uses the `Minimap2Aligner` index adapter. `relationships.py` classifies coordinate-bearing hits; `graph.py`, `decision_policy.py`, and `path_planning.py` preserve branches and defer unsupported overlap components. `merge.py` independently verifies terminal strings before appending a suffix.

The output layer is deterministic and atomic. Coordinates are zero-based, half-open internally. Every source member remains represented in provenance, including exact aliases, containments, constructed paths, and deferred nodes.
