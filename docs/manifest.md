# Manifest

The manifest is tab-separated, has one header row, and is resolved relative to its own directory.

| Column | Required | Meaning |
| --- | --- | --- |
| `sample` | yes | Unique sample identifier; output ordering is deterministic. |
| `contigs` | yes | Source FASTA path. |
| `bam` | no | Sample BAM/CRAM path; required for alignment evidence. |
| `technology` | no | Free-text technology label retained as sample metadata. |
| `assembly_graph` | no | Existing graph path; validated and retained, not used for merge authorization. |
| other columns | no | Extra values retained as metadata. |

Example:

```text
sample	contigs	bam	technology	assembly_graph
S01	assemblies/S01.fasta	alignments/S01.bam	ont	graphs/S01.gfa
S02	assemblies/S02.fasta	alignments/S02.bam	illumina	graphs/S02.gfa
```

Validate the manifest with `contigger validate --manifest samples.tsv`. Missing required columns, duplicate samples, missing files, and malformed FASTA fail before a merge starts. A missing BAM index is reported as a warning during manifest validation and becomes an error when alignment evidence is actually validated.
