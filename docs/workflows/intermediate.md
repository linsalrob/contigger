# Intermediate workflow: FASTA plus BAM/CRAM

For each sample, use a BAM/CRAM mapped against that sample's own source assembly:

```text
sample	contigs	bam	technology
S01	assemblies/S01.fasta	alignments/S01.bam	ont
S02	assemblies/S02.fasta	alignments/S02.bam	illumina
```

The alignment must be indexed and its reference names and lengths must match the FASTA exactly.

```bash
contigger validate --manifest samples.tsv
contigger validate-alignments --manifest samples.tsv
contigger merge --manifest samples.tsv --output-prefix results/contigger \
  --evidence alignments --threads 16
```

Source alignments provide sample-scoped coverage and pileup context. They do not directly contain reads spanning a junction that did not exist in the source assembly. Therefore the current merge records imperfect-junction diagnostics but does not authorize an unreviewed SNP/indel consensus; those joins remain deferred.

Create common BAMs with standard tools:

```bash
minimap2 -ax map-ont -t 16 assemblies/S01.fasta reads/S01.fastq.gz \
  | samtools sort -@ 4 -o alignments/S01.bam
samtools index alignments/S01.bam
```

For short reads, `minimap2 -ax sr` is a compatible example. Other aligners are acceptable when they produce coordinate-sorted, indexed BAM/CRAM files with matching reference names and lengths.
