# Contigger documentation

Contigger is an experimental, conservative command-line tool for reconciling redundant and overlapping contigs from related assemblies. It preserves provenance and leaves ambiguous relationships unresolved.

```text
manifest → validate → catalogue → candidates → minimap2 → relationships
         → graph decisions → safe path construction → FASTA + diagnostics
```

Start here:

- [Installation](installation.md) for Python and external tools;
- [Quick start](quickstart.md) for a FASTA-only run;
- [Beginner workflow](workflows/beginner.md) if you have assemblies only;
- [Intermediate workflow](workflows/intermediate.md) if every sample also has a mapped BAM/CRAM;
- [Advanced workflow](workflows/advanced.md) for evidence, HPC, and diagnostics.

!!! warning
    Contigger is experimental. A missed merge is preferable to a false merge. Keep original assemblies, inspect `provenance.tsv` and `ambiguous.tsv`, and do not interpret an exact duplicate as proof of biological identity.

The [inputs](inputs.md), [merge reference](usage/merge.md), [outputs](outputs.md), and [evidence guide](evidence.md) describe current behavior. Implementation and benchmark material belongs in [development](development/index.md).
