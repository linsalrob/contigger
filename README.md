# Contigger

Contigger is a conservative, provenance-aware tool for reconciling assembled metagenomic, microbial, viral, and phage contigs across samples.

> **Experimental:** Contigger is an early scaffold. It does not yet perform sequence merging or produce biological result files.

The governing principle is simple: **a missed merge is preferable to a false merge.** Identity alone never establishes that a join is valid.

## Current status

The Python 3.11+ package currently provides typed public models, transparent plain/gzip FASTA and PAF input, strict manifest validation, exact strand-aware sequence cataloguing with complete provenance, canonical positional-minimiser candidate generation, selective-alignment request planning, conservative complete-pair relationship classification, deterministic benchmark evaluation, a path-based minimap2 adapter, and a functional dry run. Graph construction, graph simplification, consensus, merging, and read analysis remain planned.

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

The experimental minimap2 path adapter supports configurable threads and `asm5`, `asm10`, or `asm20` assembly presets while recording the version and exact command. Ordinary tests and checked-in PAF classification do not require minimap2. The small fixed-seed synthetic external-tool benchmark remains available:

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

## Roadmap

1. Integrate and baseline PAF relationship classification on checked-in Pseudomonas truth. (complete)
2. Implement a stable sequence catalogue with exact and reverse-complement deduplication and complete provenance. (complete)
3. Implement canonical positional-minimiser candidates and selective-alignment planning. (complete, experimental baseline)
4. Benchmark candidate-to-alignment-to-relationship recall and add persistent minimap2 indexing/safe batching.
5. Build ambiguity-preserving containment/overlap graphs.
6. Add provenance-complete linear-path merging.
7. Add sample-aware source evidence and targeted junction remapping.

See [DESIGN.md](DESIGN.md) for assumptions, boundaries, and open questions. Contigger does **not** yet replace read-aware assembly polishing or strain-resolved assembly.
