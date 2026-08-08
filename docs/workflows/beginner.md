# Beginner workflow: FASTA only

This workflow assumes only `sample1.fasta`, `sample2.fasta`, and similar assemblies.

1. Install Python, Contigger, and minimap2 as described in [installation](../installation.md).
2. Put FASTA files under one directory and create a tab-separated manifest:

```bash
printf 'sample\tcontigs\n' > samples.tsv
for f in assemblies/*.fasta; do
  sample=$(basename "$f" .fasta)
  printf '%s\t%s\n' "$sample" "$f" >> samples.tsv
done
```

3. Validate and optionally dry-run:

```bash
contigger validate --manifest samples.tsv
contigger merge --manifest samples.tsv --output-prefix results/contigger --dry-run
```

4. Run the merge:

```bash
contigger merge --manifest samples.tsv --output-prefix results/contigger --identity 98 --threads 16
```

The useful beginner defaults are 98% identity, 1,000 bp minimum overlap, 98% containment coverage, and 50 bp end tolerance. Leave minimiser settings alone unless candidate counts or runtime motivate investigation. Sequence-only mode intentionally defers SNP/indel disagreements and ambiguous paths. Read [outputs](../outputs.md) to see what was retained and why.
