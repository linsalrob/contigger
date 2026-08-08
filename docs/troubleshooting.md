# Troubleshooting

| Symptom | Check |
| --- | --- |
| `minimap2` unavailable | Run `minimap2 --version`; install it or use `--dry-run` for input-only validation. |
| `samtools` unavailable | Install it for `validate-alignments` or `--evidence alignments`. |
| Manifest path failure | Paths are relative to the manifest, not the shell's current directory. |
| Duplicate FASTA identifier | Rename duplicates within that sample. |
| BAM reference mismatch | Re-map against the exact source FASTA; names and lengths must match. |
| BAM not indexed | Create `.bai`, `.crai`, or the supported adjacent index before validation. |
| No candidate pairs | Inspect `candidates.tsv`; sequences may be unrelated or minimiser settings too strict. |
| Too many candidates | Inspect minimiser frequency and terminal evidence before increasing resources. |
| No merges | Read `ambiguous.tsv`; high identity alone is insufficient. |
| Many deferred overlaps | This is expected for imperfect, repeat-driven, or conflicting paths. |
| Stale index rejected | Use a new `--index-dir`; indexes are content/preset/version validated. |
| Output path failure | Ensure the output parent is writable and is not an existing file. |
