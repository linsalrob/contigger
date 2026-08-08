# Contigger

Contigger conservatively combines overlapping and redundant contigs across multiple assemblies while preserving ambiguous or conflicting sequences instead of forcing unsafe joins.

It is designed for collections from metagenomes, microbial genomes, viromes, bacteriophages, and repeated or closely related assemblies. The governing rule is simple:

> **A missed merge is preferable to a false merge.**

Contigger is experimental but usable. Keep the original assemblies and inspect the provenance and ambiguity reports before using results for an important analysis.

## Why use Contigger?

Imagine three assemblies:

```text
sample_1.fasta → contigs A B C
sample_2.fasta → contigs D E F
sample_3.fasta → contigs G H I
```

Across samples, some sequences may be exact duplicates, reverse-complement duplicates, contained within longer contigs, or connected by a conflict-free terminal overlap. Others may be genuinely different, repeat-driven, strain-specific, or too ambiguous to reconcile safely. Contigger reduces redundancy when the sequence geometry supports it and leaves uncertain relationships visible for review.

It does not treat a high identity score as permission to concatenate sequences. Internal similarity is not a terminal overlap, and a BAM mapped to the original contigs cannot by itself prove a newly created junction.

## Installation

Contigger is installed from this repository; it is not currently advertised as a PyPI or Bioconda package.

```bash
git clone https://github.com/linsalrob/contigger.git
cd contigger
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

For ordinary merges, install `minimap2` in your environment. For BAM/CRAM validation or `--evidence alignments`, install `samtools` too. A convenient Mamba environment is:

```bash
mamba create -n contigger -c conda-forge -c bioconda \
  python=3.12 minimap2 samtools pip
mamba activate contigger
python -m pip install .
```

Check the installation with `contigger --version`, `contigger --help`, and, when relevant, `minimap2 --version` and `samtools --version`.

## Five-minute quick start

Put your assemblies in a directory:

```text
assemblies/
├── sample01.fasta
├── sample02.fasta
└── sample03.fasta
```

Create a tab-separated manifest named `samples.tsv`:

```text
sample	contigs
sample01	assemblies/sample01.fasta
sample02	assemblies/sample02.fasta
sample03	assemblies/sample03.fasta
```

For a directory of assemblies, Contigger can generate this manifest and discover matching GFA/BAM/CRAM sidecars:

```bash
contigger make-manifest assemblies --output samples.tsv
```

Review the generated file, then validate it before merging.

Validate it, then run the conservative sequence-only merge:

```bash
contigger validate --manifest samples.tsv
contigger merge \
  --manifest samples.tsv \
  --output-prefix results/contigger \
  --identity 98 \
  --threads 16
```

The output prefix creates:

- `results/contigger.fasta` — representatives and safely constructed sequences;
- `results/contigger.provenance.tsv` — a trace for every source contig;
- `results/contigger.relationships.tsv` — classified pair relationships;
- `results/contigger.ambiguous.tsv` — deferred components and reasons;
- `results/contigger.gfa` — graph links when `--emit-gfa` is requested;
- `results/contigger.stats.json` — configuration, counts, tools, and timings;
- `results/contigger.join_support.tsv` and `results/contigger.variants.tsv` — explicit evidence diagnostics.

Use `--dry-run` to validate inputs and print the normalized plan without writing biological outputs.

## Three levels of use

| Level | You have | Recommended mode |
| --- | --- | --- |
| Beginner | FASTA assemblies | `--evidence none` (the default) |
| Intermediate | FASTA plus BAM/CRAM mapped to each source assembly | `validate-alignments`, then `--evidence alignments` |
| Advanced | Multiple assemblies, indexed BAM/CRAM, technology metadata, graphs or raw reads | Use the supported manifest/evidence validation and inspect the separate diagnostic workflows |

`--evidence alignments` validates sample-scoped BAM/CRAM inputs and records evidence diagnostics. It does not currently authorize an unreviewed SNP/indel consensus, so imperfect overlaps remain deferred. `assembly_graph` and `technology` are accepted manifest fields for validation and provenance, but they do not independently authorize a merge. Raw FASTQ files are not direct inputs to `contigger merge`; targeted-remapping benchmark commands use checked-in datasets rather than an arbitrary raw-read manifest.

Read the [beginner workflow](docs/workflows/beginner.md), [intermediate workflow](docs/workflows/intermediate.md), or [full documentation](https://contigger.readthedocs.io/).

## What is safe to expect?

With `--evidence none`, Contigger can emit exact and reverse-complement representatives, uniquely eligible containments, and unambiguous terminal overlaps whose aligned bases are identical after orientation. It will retain both sequences when there is a SNP or indel disagreement, branch, repeat ambiguity, cycle, orientation conflict, or known-forbidden edge. A 98% identity threshold does not mean that every 98%-identical pair will merge.

Every source contig remains recoverable through provenance. Treat a deferred relationship as useful scientific information, not as a failed run.

## A practical review loop

Start with the default settings and a small representative collection. Use `--dry-run` to catch path and tool problems before creating results. After a real run, compare the FASTA count with the input count, read the top-level counts in `stats.json`, and inspect `ambiguous.tsv` before deciding whether any threshold should change. The candidate report from `contigger candidates` is useful when runtime is unexpectedly high; it tells you which pairs reached alignment, not which pairs are biologically related.

For every sequence that was removed from the representative FASTA, locate its row in `provenance.tsv`. Exact duplicates should be labeled as catalogue identities, containments should identify their surviving container, and deferred or ambiguous sequences should still have an output representative. For a constructed path, check the ordered source members, orientations, overlap coordinates, and decision reason. If a relationship matters to a biological conclusion, retain the original contigs and independently inspect the source assemblies and reads.

## Choosing an evidence mode

Use `--evidence none` when you want a reproducible sequence-only reduction or do not have mapped reads. It requires no BAM/CRAM and is the simplest starting point. Use `--evidence alignments` only when each sample has a coordinate-sorted, indexed BAM/CRAM mapped to its own source FASTA. This mode validates those references and records samtools provenance, but it remains conservative about imperfect overlaps. It is not a read-polishing command and does not pool samples into a majority consensus.

The `--minimap2-preset` option (`asm5`, `asm10`, or `asm20`) is an alignment sensitivity choice, not a biological policy. The default `asm20` is retained for current benchmark compatibility. On large collections, set `--index-dir` to a fast scratch location and monitor candidate counts, index reuse, and stage timings in `stats.json`. Do not delete the index until the run has been archived and reviewed.

## Citation, contribution, and license

A formal citation will be added when available. Contributions should follow [CONTRIBUTING.md](CONTRIBUTING.md); developer architecture and benchmark notes are in [the development documentation](docs/development/). The repository is distributed under the license in [LICENSE](LICENSE).
