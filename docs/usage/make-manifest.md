# `contigger make-manifest`

Generate a valid tab-separated manifest by scanning a directory for FASTA files and matching sidecars:

```bash
contigger make-manifest assemblies --output samples.tsv
```

The command searches recursively by default. It recognizes `.fa`, `.fasta`, `.fna`, `.fas`, and `.fsa`, including `.gz` forms. The filename with that suffix removed becomes the sample name. For each assembly it looks for matching `.gfa`/`.gfa.gz` and `.bam`/`.cram` files. It also reports whether an adjacent BAI/CRAI index was found; the index is not a manifest column.

Use `--no-recursive` to scan only the supplied directory. Duplicate sample names and multiple matching sidecars fail instead of choosing arbitrarily. Paths are written relative to the output manifest when possible, otherwise as absolute paths, so the resulting file can be passed directly to `contigger validate`.

```bash
contigger make-manifest assemblies -o samples.tsv
contigger validate --manifest samples.tsv
contigger merge --manifest samples.tsv --output-prefix results/contigger
```

The equivalent repository helper is `python scripts/make_manifest.py`; installed users should prefer the `contigger make-manifest` command.
