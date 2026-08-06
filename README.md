# Contigger

Contigger is a conservative, provenance-aware tool for reconciling assembled metagenomic, microbial, viral, and phage contigs across samples.

> **Experimental:** Contigger is an early scaffold. It does not yet perform sequence merging or produce biological result files.

The governing principle is simple: **a missed merge is preferable to a false merge.** Identity alone never establishes that a join is valid.

## Current status

The Python 3.11+ package currently provides typed public models, transparent plain/gzip FASTA and PAF input, strict manifest validation, exact strand-aware sequence cataloguing with complete provenance, canonical positional-minimiser candidate generation, pair-safe indexed alignment batching, conservative complete-pair relationship classification, deterministic benchmark evaluation, typed ambiguity-preserving relationship graph construction, conservative graph decision eligibility, provenance-complete metadata-only linear-path planning, sample-scoped source BAM/CRAM validation and pileups, targeted provisional-junction read extraction/remapping, a minimap2 adapter with validated persistent target indexes, and a functional dry run. Graph simplification, sequence construction, consensus, and merging remain planned.

## Development installation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

Validate inputs:

```bash
contigger --help
contigger validate --manifest samples.tsv
contigger merge --manifest samples.tsv --output-prefix results/contigger --dry-run
```

A real merge fails clearly because biological merging is not implemented:

```bash
contigger merge --manifest samples.tsv --output-prefix results/contigger
```

## Experimental PAF diagnostics

Classify all ordered query-target groups in a minimap2 PAF file:

```bash
contigger classify-paf --paf alignments.paf --output relationships.tsv \
  --identity 98 --min-overlap 1000 --min-containment 500 \
  --containment-coverage 98 --end-tolerance 50
```

PAF and TSV coordinates are zero-based and half-open. Blank PAF lines are ignored; malformed fields, coordinates, orientations, or tags fail with the physical source line. Mapping quality 255 is retained as valid PAF data. The TSV is diagnostic and experimental and never claims that contigs were merged.

Every distinct primary, secondary, supplementary/inversion-labelled alignment is classified geometrically. Exact duplicate records are ignored. Equivalent hits can collapse only when topology, orientation, and coordinates agree within end tolerance. Conflicting accepted hits—including repeat placements—become `AMBIGUOUS_OVERLAP`; alignment score and primary status do not elect a winner. Rejected hits and reasons remain counted.

## Manifest

The manifest is tab-separated. Only `sample` and `contigs` are currently required; relative paths are resolved from the manifest directory.

```text
sample	contigs	bam	technology	assembly_graph
S01	S01.contigs.fasta	S01.sorted.bam	illumina	S01.gfa
S02	S02.contigs.fasta	S02.sorted.bam	ont	S02.gfa
```

Unknown columns are retained as sample metadata. Optional paths, when supplied, must exist.

## Planned outputs

The first implemented outputs will be `contigger.fasta`, `contigger.provenance.tsv`, `contigger.relationships.tsv`, `contigger.ambiguous.tsv`, `contigger.gfa`, and `contigger.stats.json`. Later milestones plan `contigger.variants.tsv`, `contigger.join_support.tsv`, `contigger.consensus.vcf`, and `contigger.low_confidence.bed`. A dry run writes none of these.

## External tools and synthetic benchmark

The experimental minimap2 adapter supports configurable threads and `asm5`, `asm10`, or `asm20` assembly presets while recording the version and exact command. Persistent target indexes are built with the selected preset and carry deterministic metadata containing that preset, the minimap2 version, target identifiers, lengths, and a sequence checksum; incomplete or mismatched indexes are rejected. Selective batches contain multiple approved queries but exactly one target, preventing accidental all-v-all expansion. Ordinary tests and checked-in PAF classification do not require minimap2. The small fixed-seed synthetic external-tool benchmark remains available:

When minimap2 is installed, compare `asm5` and `asm20` against fixed-seed synthetic truth:

```bash
python benchmarks/evaluate_minimap2.py
python benchmarks/evaluate_minimap2.py --json benchmark.json
pytest -m integration
```

The report puts false merges first, then correct, missed, ambiguous, candidate-record, and timing counts. The fixtures cover orientation, containment, terminal overlaps, substitutions, an indel, internal similarity, repeats, incompatible placements, low complexity, unrelated sequence, and classifier boundaries. This small benchmark is experimental and does not establish production sensitivity.

## Checked-in Pseudomonas benchmark

`test_data` version 1.0.0 contains 90 derived contigs across three samples and 74 construction-derived ordered-pair truth rows. It includes checked-in `asm5` and `asm20` PAFs, but not the full source reads, assembly FASTA, or GFA. Evaluate either PAF without minimap2:

```bash
contigger benchmark --dataset test_data \
  --paf test_data/alignments/all_vs_all.asm20.paf.gz \
  --output-json benchmark.json --output-tsv benchmark.tsv
```

A **false merge** is a merge-like pair result where pairwise truth forbids merging. A **missed relationship** is an absent or `NO_RELATIONSHIP` result for an unambiguous valid truth pair. Unexpected PAF pairs are separate because the construction table is sparse. Multiple incompatible hits within one ordered pair are pair-level ambiguity; ambiguity across different targets requires graph/component context and is explicitly deferred.

At the default thresholds, `asm5` has 58 correct classifications, 2 false merges, and 4 missed relationships; `asm20` has 60 correct, 2 false merges, and 2 missed relationships. Both false merges are the two ordered directions of the 51 bp end-tolerance case: minimap2 extends through identical flanking sequence, so the PAF geometry looks terminal. `asm20` is more sensitive on this dataset, but the production default remains configurable and unchanged. Broader microbial, viral, and phage truth sets are required before changing it. The deterministic baseline is in `benchmarks/pseudomonas_baseline.json`.

This command scores classifier output only. It does not construct a graph or claim that any contigs were merged.

## Exact catalogue and candidate planning

Create deterministic canonical sequences while retaining every source contig and explicit strand in provenance:

```bash
contigger catalogue --manifest samples.tsv \
  --output-fasta catalogue.fasta \
  --output-provenance catalogue.provenance.tsv
```

Only byte-exact forward or reverse-complement sequences collapse. Catalogue identifiers are derived from the SHA-256 digest of the lexicographically canonical strand; this is exact deduplication, not biological merging.

Generate typed positional-minimiser evidence for selective alignment:

```bash
contigger candidates --manifest samples.tsv --output candidates.tsv \
  --kmer-size 21 --window-size 10 --min-shared-minimisers 5 \
  --max-minimiser-frequency 100 --terminal-band 1000
```

Ambiguous-base k-mers and globally frequent minimisers are excluded. Shared minimisers alone are not emitted as candidates: retained pairs must have positional evidence at compatible sequence ends or across both ends of a possible contained sequence. Multiple orientation/topology signals remain explicit. Candidate rows request further alignment; they are not relationships and cannot authorize a merge.

On Pseudomonas benchmark 1.0.0, exact/RC deduplication reduces 90 sources to 84 canonical sequences. The documented candidate settings retain all 23 valid non-exact truth case groups in 61 candidates out of 3,486 possible canonical pairs. This is a recall/integration baseline, not evidence that all 61 pairs are mergeable; forbidden boundary and ambiguous repeat pairs deliberately proceed to alignment and conservative classification.

Evaluate the complete deterministic pathway using the checked-in PAF as alignment observations (minimap2 is not invoked):

```bash
contigger benchmark-pipeline --dataset test_data \
  --paf test_data/alignments/all_vs_all.asm20.paf.gz \
  --output-json pipeline-benchmark.json
```

Both presets recover all 12 exact/RC truth rows and all 23 valid non-exact case groups (40 ordered truth rows); candidate generation adds zero misses. At relationship classification, both retain the two known 51 bp end-tolerance false merges, while `asm5` misses four relationships and `asm20` misses two. The pair-stage report defers four graph-level ambiguity groups. The checked-in result is `benchmarks/pseudomonas_pipeline_baseline.json`. This is a staged recall and safety benchmark, not a graph result or evidence that any contigs were merged.

## Ambiguity-preserving relationship graphs

`build_relationship_graph()` consumes complete `PairRelationship` decisions and returns stable nodes, structurally separate containment and terminal-overlap edges, explicit ambiguous edges, and deterministically ordered components. Consistent reciprocal PAF decisions collapse to one edge; inconsistent reciprocals remain ambiguous. Components are marked ambiguous for competing use of one oriented terminal, multiple possible containers, pair-level ambiguity, cycles, or contradictory orientations. A degree-two linear path is retained without being declared merged.

Exact matches are rejected at this boundary because exact and reverse-complement deduplication belongs in the catalogue stage. `NO_RELATIONSHIP` decisions do not create edges, although supplied isolated nodes remain in the graph. Graph construction performs no containment removal, edge selection, path simplification, sequence joining, or merge authorization.

The source-identifier diagnostic regression in `benchmarks/pseudomonas_graph_baseline.json` has 90 nodes and three containment edges for both PAFs. `asm5` has 50 overlap edges in 58 components; `asm20` has 51 in 57. In both, the four deferred repeat/branch truth groups occur in one preserved ambiguous component. The known forbidden 51 bp end-tolerance pair also remains an ordinary edge, demonstrating why graph presence alone cannot authorize merging.

`evaluate_graph_decisions()` is the next conservative boundary. It marks a unique containment in an unambiguous component eligible for later provenance-aware disposition, but does not remove the contained node. An overlap-only component is eligible for later path planning only when it is unambiguous and every proposed junction edge has explicit support supplied by a future evidence stage. Containment-mixed or ambiguous components stay deferred. Eligibility is not merge authorization and produces no biological output.

With no junction evidence supplied, `benchmarks/pseudomonas_decision_policy_baseline.json` records three eligible containment dispositions for each preset and zero eligible overlap components. All 50 `asm5` and 51 `asm20` overlap edges remain deferred, including both the repeat-connected ambiguous component and the known forbidden 51 bp boundary edge.

`plan_linear_paths()` accepts a complete catalogue and matching relationship graph, invokes the conservative decision policy, and returns canonical path metadata only for eligible overlap components. Each node carries explicit path orientation and every exact/RC catalogue source member with its path-relative orientation. Reverse-complement-equivalent traversals collapse to one deterministic representation. The planner rejects graph/catalogue mismatches and never trims, joins, or writes sequence.

The checked-in source-ID diagnostic path baseline supplies no junction evidence, so both PAF presets produce zero paths. All 19 `asm5` and 20 `asm20` overlap components remain deferred.

## Source alignment evidence

Validate supplied BAM/CRAM files against their sample FASTA references using installed `samtools`:

```bash
contigger validate-alignments --manifest samples.tsv
```

`BamEvidenceProvider` requires an adjacent BAI/CRAI, checks file integrity and index readability, and requires exact `@SQ` reference names and lengths. Pileups accept zero-based half-open source intervals and return sample-labelled allele counts, depth, and mean base/mapping qualities. These observations describe existing source contigs only; they cannot validate a newly constructed junction.

`TargetedJunctionRemapper` tests an explicitly supplied provisional reference rather than constructing one. It uses `samtools view` to nominate primary read names near two named source-contig ends, recovers all primary records (including mates) sharing those names, name-collates and converts the bounded subset to FASTQ, and remaps it with minimap2. A distinct read counts as spanning only when its primary alignment crosses the supplied zero-based junction coordinate by the configured number of reference bases on both sides. Results retain sample identity, exact tool versions and commands, selected/remapped/spanning read names, and diagnostics.

This report is evidence, not merge authorization. It does not infer trimming, construct joined sequence, call consensus, pool samples, or automatically supply graph decision-policy support. Minimum read counts, mapping-quality rules, technology-specific presets, duplicate handling, and adequate flank lengths require reviewed benchmarks before junction observations can affect a biological decision.

## Roadmap

1. Integrate and baseline PAF relationship classification on checked-in Pseudomonas truth. (complete)
2. Implement a stable sequence catalogue with exact and reverse-complement deduplication and complete provenance. (complete)
3. Implement canonical positional-minimiser candidates and selective-alignment planning. (complete, experimental baseline)
4. Benchmark candidate-to-alignment-to-relationship recall and add persistent minimap2 indexing/safe batching. (complete, experimental baseline)
5. Build typed, ambiguity-preserving containment/overlap graphs without simplifying or merging them. (complete, unsimplified experimental baseline)
6. Define and benchmark conservative graph decision policies for containment disposition and merge-path eligibility. (complete; no biological output)
7. Implement provenance-complete unambiguous linear-path planning without sequence merging. (complete; metadata only)
8. Validate sample-aware source BAM/CRAM references and expose source-contig evidence without junction claims. (complete)
9. Add targeted junction read extraction, remapping, and evidence reporting. (complete; evidence only)
10. Benchmark technology-specific junction-support and consensus/variation policies before connecting evidence to graph decisions.

See [DESIGN.md](DESIGN.md) for assumptions, boundaries, and open questions. Contigger does **not** yet replace read-aware assembly polishing or strain-resolved assembly.
