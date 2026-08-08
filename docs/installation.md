# Installation

## Requirements

Contigger requires Python 3.11 or newer. `minimap2` is required for a normal merge whenever candidate pairs need alignment. `samtools` is required for BAM/CRAM validation and for `--evidence alignments`. Unit tests and dry runs do not need either executable.

## Recommended Mamba installation

```bash
mamba create -n contigger -c conda-forge -c bioconda \
  python=3.12 minimap2 samtools pip
mamba activate contigger
git clone https://github.com/linsalrob/contigger.git
cd contigger
python -m pip install .
```

If you only have FASTA files, `samtools` is optional. Verify:

```bash
contigger --version
contigger --help
minimap2 --version
samtools --version  # needed for BAM/CRAM workflows
```

## Developer and documentation installs

```bash
python -m pip install -e '.[dev]'
python -m pip install -e '.[docs]'
```

The first installs pytest, Ruff, and mypy. The second installs MkDocs and the Material theme; no mapper or read tool is installed by the documentation extra.
