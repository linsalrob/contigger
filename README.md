# Contigger

Contigger is a conservative, provenance-aware tool for reconciling assembled metagenomic, microbial, viral, and phage contigs across samples.

> **Experimental:** Contigger is an early scaffold. It does not yet perform sequence merging or produce biological result files.

The governing principle is simple: **a missed merge is preferable to a false merge.** Identity alone never establishes that a join is valid.

## Current status

The Python 3.11+ package currently provides typed public models, strict FASTA and TSV-manifest validation, normalised run configuration, a pure conservative relationship classifier, deterministic output naming/provenance interfaces, minimap2 and BAM/CRAM adapter boundaries, and a functional dry run. Graph simplification, consensus, merging, read analysis, and production candidate generation remain planned.

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

## External tools

Future workflows may wrap minimap2 for assembly alignment and targeted remapping, and samtools for BAM/CRAM validation, indexing, extraction, and conversion. Neither tool is required for unit tests or input-only dry runs. `pysam` is not a mandatory dependency.

## Roadmap

1. Complete and benchmark PAF relationship classification on synthetic contigs.
2. Implement positional-minimiser candidate generation.
3. Build ambiguity-preserving containment/overlap graphs.
4. Add provenance-complete linear-path merging.
5. Add sample-aware source evidence and targeted junction remapping.

See [DESIGN.md](DESIGN.md) for assumptions, boundaries, and open questions. Contigger does **not** yet replace read-aware assembly polishing or strain-resolved assembly.
