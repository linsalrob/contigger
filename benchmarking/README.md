# Contigger deterministic benchmark generator

This directory builds a compact validation dataset from one circular assembly, its GFA, and ONT reads. It deliberately applies the rule **a missed merge is preferable to a false merge**. Complete source data are never copied into the output.

## Requirements and environment

Activate the Mamba environment containing Python 3, PyYAML, pysam, pytest, minimap2, samtools, seqkit and pigz. Inputs default to `BC04.trimmed.fastq.gz`, `consensus_assembly.fasta`, and `consensus_assembly.gfa` in the parent/workspace directory. No downloads are performed.

## Build, validation and regeneration

```bash
python benchmarking/build_benchmark.py --config benchmarking/config.yaml --output test_data
python benchmarking/validate_benchmark.py test_data
```

Use `--force` to replace an existing output, `--seed` and the input path flags to override configuration, and `--keep-work` to retain private intermediates. Otherwise `benchmarking/.work/` (the minimap2 index, source self-alignment, whole-genome BAM, and combined FASTA/PAFs) is deleted after success. Fixed seeds, sorted records, gzip timestamp zero, and construction-defined truth make regeneration deterministic.

The configuration defines inclusive identity (98%), overlap (1,000 bp), containment coverage (98%), and end-tolerance (50 bp) thresholds. Cases cover exact/RC equality, directional containment, both overlap directions, RC overlaps, SNP/indel error, identity/length/end/coverage boundaries, circular origin, unrelated/internal/chimeric/low-complexity negatives, and controlled repeat ambiguity. Synthetic repeats are labelled honestly when no natural repeat is asserted.

Every file is limited to 50 MiB; the preferred archive size is 100 MiB and hard maximum is 150 MiB. Targeted read retention can be reduced deterministically with `reads.maximum_reads_per_locus` if necessary.

To use in Contigger, copy `benchmarking/` and `test_data/` (and optionally `test_data.tar.gz`) into the repository root. The source FASTQ, assembly, GFA, whole-genome BAM, and indexes are intentionally excluded.
