# Validation commands

Validate a manifest and all referenced FASTA/optional files:

```bash
contigger validate --manifest samples.tsv
```

Validate sample BAM/CRAM indexes, integrity, and exact FASTA references:

```bash
contigger validate-alignments --manifest samples.tsv
```

The second command requires `samtools` and at least one manifest BAM/CRAM. Paths are resolved relative to the manifest. Warnings such as a missing adjacent alignment index are reported during ordinary validation; alignment validation fails until the index exists.
