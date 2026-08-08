# Reproducibility

Keep the input manifest, exact Contigger version/commit, command line, Python environment, minimap2 and samtools versions, and the complete output prefix together. `stats.json` records normalized configuration, tool versions and commands, index metrics, counts, and stage timings. `provenance.tsv` records source membership and coordinates; `relationships.tsv` and `ambiguous.tsv` preserve decision diagnostics.

Do not overwrite original assemblies or silently reuse an index from a different sequence set, preset, or minimap2 version. Contigger validates index metadata and rejects mismatches. For published analyses, archive checksums of input FASTA/BAM/CRAM files and the repository commit.
