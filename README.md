# Contigger

Contigger is a conservative, provenance-aware tool for reconciling assembled metagenomic, microbial, viral, and phage contigs across samples.

> **Experimental:** Contigger is an early scaffold. It does not yet perform sequence merging or produce biological result files.

The governing principle is simple: **a missed merge is preferable to a false merge.** Identity alone never establishes that a join is valid.

## Current status

The Python 3.11+ package currently provides typed public models, strict FASTA and TSV-manifest validation, normalised run configuration, strict PAF stream parsing, conservative complete-pair relationship classification, a path-based minimap2 adapter, deterministic output/provenance interfaces, and a functional dry run. Graph construction, graph simplification, consensus, merging, read analysis, and production candidate generation remain planned.

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

The experimental minimap2 path adapter supports configurable threads and `asm5`, `asm10`, or `asm20` assembly presets while recording the version and exact command. Ordinary tests and PAF classification do not require minimap2. The existing `asm20` default remains unchanged: minimap2 was unavailable in the development environment, so there is not yet local evidence sufficient to replace it despite the conceptual 98% identity threshold.

When minimap2 is installed, compare `asm5` and `asm20` against fixed-seed synthetic truth:

```bash
python benchmarks/evaluate_minimap2.py
python benchmarks/evaluate_minimap2.py --json benchmark.json
pytest -m integration
```

The report puts false merges first, then correct, missed, ambiguous, candidate-record, and timing counts. The fixtures cover orientation, containment, terminal overlaps, substitutions, an indel, internal similarity, repeats, incompatible placements, low complexity, unrelated sequence, and classifier boundaries. This small benchmark is experimental and does not establish production sensitivity.

## Roadmap

1. Complete and benchmark PAF relationship classification on synthetic contigs.
2. Implement positional-minimiser candidate generation.
3. Build ambiguity-preserving containment/overlap graphs.
4. Add provenance-complete linear-path merging.
5. Add sample-aware source evidence and targeted junction remapping.

See [DESIGN.md](DESIGN.md) for assumptions, boundaries, and open questions. Contigger does **not** yet replace read-aware assembly polishing or strain-resolved assembly.
