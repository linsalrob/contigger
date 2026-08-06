# Contigger initial design

## 1. Problem statement

Assemblies from related samples often contain exact duplicates, contained sequences, compatible terminal overlaps, strain-specific alternatives, repeats, and errors. Contigger will reconcile them without converting similarity into unjustified joins. Its priority is recoverability and avoidance of false merges, not maximal compression.

## 2. Terminology

A **contig** is an input assembled sequence. A **candidate** is a pair selected by positional seed evidence. An **alignment hit** is a coordinate-bearing observation, not a decision. **Containment** means one sequence is covered end-to-end within another. A **terminal overlap** connects compatible ends. A **relationship** is the classified interpretation of an alignment. A **junction** is newly constructed sequence adjacency. **Provenance** maps every source interval and disposition to an output or ambiguity record.

## 3. Scope of the current implementation

The current milestone supplies typed models, transparent plain/gzip FASTA and PAF input, manifest validation, stable exact strand-aware sequence cataloguing, canonical positional-minimiser candidates, selective per-pair alignment requests, strict streaming PAF parsing, conservative classification of complete ordered query-target hit groups, deterministic evaluation against checked-in Pseudomonas truth, and unsimplified ambiguity-preserving relationship graph construction. Catalogue, candidate, and graph objects are evidence/planning artifacts; no current command simplifies a graph or merges sequence.

## 4. Explicit non-goals

This scaffold does not implement graph simplification, containment removal, merge-path authorization, path merging, consensus construction, targeted read remapping, variant calling, confidence scoring, parallel or distributed execution, or Rust bindings. It must never imply that these operations succeeded.

## 5. Input model

Each `SampleInput` has a unique sample identifier and contig FASTA, with optional coordinate-sorted indexed BAM/CRAM, technology, assembly graph, and extra metadata. FASTA IDs need only be unique within a sample; internal identities include sample context. Optional inputs are never silently ignored when supplied.

## 6. Manifest format

The tab-separated manifest requires `sample` and `contigs`. Optional columns are `bam`, `technology`, and `assembly_graph`; unknown columns are preserved as strings in metadata. Fields are trimmed, blank required fields and duplicate samples are errors, paths resolve relative to the manifest, and diagnostics contain line numbers.

## 7. Internal coordinate and orientation conventions

All internal intervals are zero-based and half-open: `[start, end)`. PAF already follows this convention; formats that do not must convert at their boundary. Coordinates remain in the named sequence's forward coordinate system. Orientation is always an explicit `+` or `-` and is never inferred from coordinate order. For reverse hits, target terminal semantics are swapped explicitly during topology checks.

## 8. Processing pipeline

The planned stages are validation; stable sequence identification; exact deduplication; positional-minimiser candidate generation; selective alignment; relationship classification; ambiguity-preserving graph construction; confident containment removal; conservative linear-path merging; optional sample-aware evidence analysis; targeted remapping for new junctions; and deterministic outputs with run provenance. Every stage consumes and produces typed records.

## 9. Candidate generation

Initial conceptual defaults, all requiring benchmarking, are:

```text
identity threshold:              98%
minimum terminal overlap:        1,000 bp
minimum containment span:        500 bp
containment coverage:            98%
end tolerance:                   50 bp
canonical k-mer size:            21
minimiser window:                10
minimum shared minimisers:       5
maximum minimiser frequency:     100
```

General k-mer composition vectors are not preferred because they discard position and can conflate shared composition with overlap. Canonical positional minimisers are the intended first heuristic. Candidate records should retain shared-seed counts, relative positions, and orientation signals. Frequent minimisers must be bounded to limit repeats and quadratic pair expansion.

Exact catalogue sequences use the lexicographically smaller forward/reverse-complement strand and a full SHA-256 content identifier. Every source maps to the full canonical interval with explicit orientation and representative/duplicate disposition. No non-exact sequence is collapsed.

Minimiser generation skips k-mers containing ambiguity symbols, hashes canonical k-mers deterministically, retains tied window minima, and globally suppresses observations above the configured frequency. Candidate evidence retains positions, every supported orientation, and every compatible terminal or possible-containment topology. End evidence may be contributed by different seeds along the same pair; requiring one seed to touch both ends would miss long overlaps. These candidates are intentionally over-inclusive alignment requests, not merge decisions.

## 10. Alignment abstraction and PAF handling

An `Aligner` can build/reuse an index, align typed sequences, return `AlignmentHit` objects, and report its name, version, and exact last command. The minimap2 adapter can safely align target and query FASTA paths with an argument array and parse the captured PAF stream. It records the executable version and exact alignment command; failures retain stderr. The backend-neutral hit model retains optional alignment role and the `AS`, `cm`, `s1`, and `s2` observations without requiring them.

Selective-alignment planning resolves only canonical candidate pairs into validated requests. The simple executor retains one-query/one-target calls. The indexed executor groups approved queries by exactly one target, builds or reuses that target index, and rejects returned identifiers outside the request set. It never sends multiple targets in a batch, so batching cannot silently expand into all-v-all alignment.

Minimap2 indexes are built with the selected assembly preset, written through a temporary same-filesystem path, and identified by deterministic sidecar metadata containing format version, preset, minimap2 version, ordered target identifiers and lengths, and a checksum over identifiers and sequences. Reuse requires an exact metadata match; incomplete, corrupt, or stale indexes fail rather than being silently trusted or overwritten. Index filenames are hashes of target identifiers, so identifiers cannot alter the index directory layout.

PAF input is processed line by line. Only blank lines are ignored. Required fields, coordinates, orientation, mapping quality (including the valid unknown value 255), and optional-tag syntax are validated; errors contain physical line numbers. Coordinates remain zero-based and half-open.

## 11. Relationship classification

`RelationshipType` defines `EXACT_MATCH`, `QUERY_CONTAINED_IN_TARGET`, `TARGET_CONTAINED_IN_QUERY`, `QUERY_SUFFIX_TO_TARGET_PREFIX`, `TARGET_SUFFIX_TO_QUERY_PREFIX`, `AMBIGUOUS_OVERLAP`, and `NO_RELATIONSHIP`. Classification considers identity, aligned length, query and target coverage, distance to every query and target end, explicit orientation, alignment topology, and every distinct hit for the ordered pair.

Single-hit classification remains the geometric primitive. At pair level, exact duplicate records are removed deterministically and rejected-hit diagnostics are retained. Primary, secondary, and inversion-labelled hits receive the same geometric scrutiny: labels and scores never select a winner. Accepted hits collapse only when relationship type, terminal topology, orientation, and all placement coordinates agree within the configured end tolerance. Incompatible types, topologies, orientations, or materially different placements produce `AMBIGUOUS_OVERLAP`. This deliberately turns unresolved repeat evidence into ambiguity rather than a score-selected relationship.

## 12. Containment criteria

Containment is evaluated separately from overlap. The contained sequence must be covered at both ends within tolerance, meet minimum containment span and coverage, and meet identity. The containing sequence need not terminate at the alignment. Near-equal sequences require deterministic representative selection and complete provenance. Threshold defaults require benchmarking.

## 13. Terminal-overlap criteria

A merge candidate must meet identity and minimum aligned span and connect exactly one compatible pair of oriented terminals within end tolerance. Internal local similarity is `NO_RELATIONSHIP`, regardless of identity. Simultaneously compatible or inconsistent topologies become ambiguous rather than being selected greedily. Sequence identity alone is never merge evidence.

## 14. Graph representation

The relationship graph has stable sequence nodes and typed relationship edges retaining orientation, representative zero-based half-open coordinates, accepted/rejected hit counts, and diagnostics. Containment, terminal overlap, and ambiguous evidence are separate edge collections. `NO_RELATIONSHIP` creates no edge, while supplied isolated nodes remain. Exact matches are rejected because they must already have been resolved by the catalogue.

Complete ordered-pair decisions are the graph boundary because bare relationships do not retain representative coordinates or competing-hit evidence. Reciprocal decisions collapse only when relationship topology, explicit orientation, identity, and swapped coordinates agree. For reverse-complement overlaps, reversing the ordered pair preserves the query-relative topology name because the same physical terminals remain involved. Reciprocal disagreement becomes an ambiguous edge; scores never choose a winner.

Components are deterministic connected components over every retained edge. They expose ambiguity when they contain pairwise ambiguous evidence, multiple possible containers for one sequence, competing overlap edges using the same physical oriented terminal, an overlap cycle, or inconsistent relative orientations. Ordinary degree two is not itself ambiguous: a linear path uses a prefix and suffix once each. Graph construction only records this structure and never removes a node or selects a path.

## 15. Conservative graph simplification

Confident contained nodes may be removed from representative output only after provenance is attached. Only unambiguous degree-compatible linear overlap paths may be proposed for merging. Branches, orientation conflicts, competing edges, repeat signals, cycles without a justified rule, and conflicting evidence are preserved. No greedy best-score collapse is permitted.

## 16. Sequence conflict handling

Evidence collection and decision policy are separate interfaces. Evidence modes are `none`, `alignments`, and `reads`. Proposed policies are `representative`, `majority`, `quality-weighted`, `sample-aware`, `ambiguous`, and `reject`. These names establish configuration boundaries only; unsupported biological decisions are not implemented. The conservative default is `reject`.

## 17. BAM/CRAM evidence model

Three claims must remain distinct: sequence-overlap evidence, source-contig pileup evidence, and support for a newly constructed join. Existing BAM/CRAM records can provide depth, allele counts, base qualities, mapping qualities, strand support, clipping, paired-read metadata, and sample identity within source contigs. Providers are sample-scoped and must validate reference names and lengths against the source FASTA.

## 18. Sample-aware variation

Evidence is retained by sample. For example:

```text
sample A: allele A strongly supported
sample B: allele G strongly supported
```

This must not become a pooled majority call without an explicit policy; it may be genuine strain or population variation. Outputs must permit both observations to remain recoverable.

## 19. Targeted remapping for new junctions

The merged sequence did not exist when the original BAM was produced, so that BAM cannot directly contain an alignment spanning its new junction. A later workflow should (1) identify reads near relevant contig ends, (2) recover mates and soft-clipped reads where possible, (3) construct a provisional merged reference, (4) remap only the relevant subset, and (5) assess junction-spanning support. Minimap2 may perform remapping and samtools may extract/convert reads; neither implies a junction decision by itself.

## 20. Provenance model

Every retained, contained, merged, ambiguous, or rejected source contig remains represented. Rows record output ID, sample, source ID, relationship, explicit orientation, zero-based half-open source and output intervals, identity where applicable, disposition, and reason. Stable output IDs must not depend on traversal timing. Writers use fixed columns and sorted rows.

## 21. Output formats

The experimental `classify-paf` command emits a deterministic diagnostic relationships TSV with zero-based half-open coordinates, accepted/rejected counts, and reasons. It is not a merge result or a stable production schema. Planned biological outputs remain `contigger.fasta`, `contigger.provenance.tsv`, `contigger.ambiguous.tsv`, `contigger.gfa`, and `contigger.stats.json`; they are not emitted by the current implementation. A dry run writes no biological outputs.

## 22. Determinism and reproducibility

Samples, candidates, relationships, graph nodes/edges, and output rows use documented stable sort keys. Configuration, exact external commands, executable versions, and input identities belong in run statistics. Future randomness must accept and record a seed. Hash-table and traversal order cannot affect output identifiers or decisions.

### 22.1 Benchmark truth and ambiguity scope

Benchmark truth parsing is separate from production `Relationship` models. It validates typed TSV fields with physical line numbers, rejects duplicate ordered pairs, and never infers truth from names or PAF output. The table is sparse, so observed pairs absent from truth are unexpected rather than negative truth.

Pair-level ambiguity means non-equivalent alignments compete within one ordered query-target group; `classify_pair()` preserves it. Graph-level ambiguity means one contig has plausible relationships to multiple targets. A pair classifier cannot observe that context, so the evaluator defers the four graph groups (`incompatible_placements`, `opposite_orientations`, `repeat_ambiguity`, and `terminal_repeat`).

### 22.2 Pseudomonas 1.0.0 baseline

At the documented defaults, `asm5` produces 58 correct classifications, 2 false merges, and 4 missed relationships. `asm20` produces 60 correct, 2 false merges, and 2 missed relationships. Both false merges are ordered directions of `end_tolerance_51`: minimap2 extends through 51 identical flanking bases, hiding the construction boundary from the pair classifier. No classifier rule changed because the PAF contains insufficient evidence for a general safe distinction. `asm20` is more sensitive here, but the default remains configurable; broader microbial, viral, and phage truth sets are needed before changing it.

### 22.3 Candidate baseline

The 90 source contigs collapse into 84 exact canonical sequences. With `k=21`, window 10, five shared minimisers, maximum global frequency 100, and a 1,000 bp terminal band, all 23 valid non-exact truth case groups reach selective alignment in 61 candidates out of 3,486 possible pairs. Candidate inclusion is not a positive relationship call: threshold-negative and ambiguity cases are expected to remain until alignment/classification. Internal-only and low-complexity negatives have permanent focused regressions.

### 22.4 Catalogue-to-relationship baseline

The stage-aware evaluator uses the checked-in complete PAFs as deterministic alignment observations and does not invoke minimap2. Both presets recover all 12 exact/RC truth rows through the catalogue and route all 23 valid non-exact case groups (40 ordered rows) through candidate selection, so candidate-stage missed relationships are zero. The eligible relationship results retain the established classifier baseline: two false merges for each preset, four missed relationships for `asm5`, two for `asm20`, and four graph-level ambiguity groups deferred. No classifier rule or threshold changed. This separation identifies which stage loses recall without treating a candidate as a relationship or merge.

### 22.5 Unsimplified graph baseline

The checked-in PAF graph regression uses source identifiers as a diagnostic proxy; production candidate alignments use catalogue identifiers. Exact/self decisions are excluded because catalogue construction precedes the graph. Both presets retain 90 nodes and three containment edges. `asm5` produces 50 overlap edges and 58 components; `asm20` produces 51 and 57. Each has one repeat-connected ambiguous component containing all four previously deferred graph ambiguity groups. The known forbidden `end_tolerance_51` pair is also retained as one overlap edge. This is intentional evidence that graph membership is not merge permission and that a later decision policy must remain benchmarked against forbidden edges.

## 23. Performance and scalability

Streaming FASTA parsing limits parser overhead, while positional minimisers and frequency filtering should avoid all-v-all alignment. Index reuse, batching, and bounded candidate sets precede parallelism. Profiling and representative benchmarks must justify optimisation. Stable typed boundaries allow hot paths to move to Rust without changing public records or the CLI.

## 24. Error handling

Malformed syntax, unreadable supplied files, missing references, duplicate identifiers, invalid coordinates, inconsistent lengths, unsupported policies, and failed tools are explicit domain errors. External failures include stderr and exact arguments. Warnings cover non-fatal conditions such as a supplied BAM lacking an adjacent index. Placeholder methods raise `FeatureNotImplementedError`.

## 25. Testing strategy

Unit tests use synthetic FASTA and manually built alignment records and require no external tools. They cover both orientations, containment directions, terminal geometries, internal similarity, thresholds, zero spans, ambiguity, malformed input, determinism, and CLI status codes. External-tool tests are marked integration tests. False-join regressions are permanent; random tests use fixed seeds and temporary files use pytest fixtures.

## 26. Staged roadmap

1. Integrate and baseline relationship classification using synthetic and checked-in Pseudomonas truth. (complete)
2. Implement a stable sequence catalogue with exact and reverse-complement deduplication and complete provenance. (complete)
3. Implement and benchmark canonical positional-minimiser candidates and selective requests. (complete, experimental baseline)
4. Benchmark executable candidate-to-alignment-to-relationship recall and add persistent indexing/safe batching. (complete, experimental baseline)
5. Add typed, ambiguity-preserving graph construction and containment handling, without graph simplification or sequence merging. (complete, unsimplified experimental baseline)
6. Define and benchmark conservative containment-disposition and merge-path eligibility policies, preserving every branch, conflict, cycle, and forbidden-edge regression.
7. Implement provenance-complete unambiguous linear-path planning before any sequence merging.
8. Validate source BAM/CRAM references and expose sample-aware pileups.
9. Add targeted read extraction/remapping and junction-support reporting.
10. Evaluate consensus and variation policies only with reviewed benchmarks.

## 27. Open design questions

- Which minimap2 preset best balances sensitivity and repeat safety on broader reviewed truth sets?
- Are containment thresholds length- or technology-dependent?
- How should circular contigs be represented without unsafe linearisation?
- Which graph patterns can be simplified safely in the presence of reverse complements?
- What constitutes adequate junction support for each sequencing technology?
- How should read names, mates, and sample metadata be retained without excessive storage?
- Which representative-selection policy is stable, explainable, and biologically neutral?
- Which benchmark truth sets adequately measure false merges across microbes, viruses, and phages?
