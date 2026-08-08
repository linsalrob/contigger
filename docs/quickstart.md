# Quick start

For three FASTA assemblies, create a tab-separated `samples.tsv`:

```text
sample	contigs
S01	assemblies/S01.fasta
S02	assemblies/S02.fasta
S03	assemblies/S03.fasta
```

Paths are resolved relative to the manifest. Run:

```bash
contigger validate --manifest samples.tsv
contigger merge --manifest samples.tsv --output-prefix results/contigger --dry-run
contigger merge --manifest samples.tsv --output-prefix results/contigger --identity 98 --threads 16
```

The dry run prints normalized settings and available external tools without creating outputs. The real run writes FASTA, provenance, relationships, ambiguity, stats, and diagnostic files. Start interpretation with [outputs](outputs.md) and [interpretation](interpretation.md).

The default is deliberately conservative: exact/RC identities, safe containment, and exact conflict-free terminal overlaps may be emitted; imperfect overlaps are deferred.
