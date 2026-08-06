# Contributing to Contigger

Contigger welcomes narrow, tested contributions that preserve its conservative biological policy.

## Local setup and verification

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
pytest
ruff check .
ruff format --check .
mypy src
```

Use `ruff format .` to format changes. Add or update tests and documentation in the same change. Unit tests must stay small, deterministic, and independent of optional command-line tools; mark external-tool integration tests explicitly.

## External tool integrations

Implement the typed protocol first, isolate the adapter, discover the executable explicitly, capture its version and exact argument array, preserve stderr in failures, and never invoke a shell. Provide mock-driven unit tests plus a marked integration test. Keep optional dependencies isolated.

## Changing merge criteria

Open a proposal that identifies the biological assumption, describes expected false-positive and false-negative effects, cites or supplies benchmark data, and explains effects on provenance and ambiguity. Changes affecting biological decisions require rationale and regression tests, including permanent tests for every known false join. Do not make a threshold more aggressive without updating `DESIGN.md` and demonstrating that the non-negotiable invariants still hold.

For any merge-criterion or relationship-classifier change, run both checked-in evaluations and explain every difference from `benchmarks/pseudomonas_baseline.json`:

```bash
contigger benchmark --dataset test_data --paf test_data/alignments/all_vs_all.asm5.paf.gz
contigger benchmark --dataset test_data --paf test_data/alignments/all_vs_all.asm20.paf.gz
```

Catalogue and candidate-generation changes must also preserve exact/RC provenance tests and the checked-in candidate baseline. Candidates are permitted to be over-inclusive, but any loss of a valid truth case requires biological review before changing `benchmarks/pseudomonas_candidates_baseline.json`.

Changes to catalogue, candidate, alignment planning, or relationship stages must run both staged evaluations and explain any difference from `benchmarks/pseudomonas_pipeline_baseline.json`:

```bash
contigger benchmark-pipeline --dataset test_data --paf test_data/alignments/all_vs_all.asm5.paf.gz
contigger benchmark-pipeline --dataset test_data --paf test_data/alignments/all_vs_all.asm20.paf.gz
```

Graph changes must run `pytest tests/test_graph.py tests/test_graph_benchmark.py` and explain every change from `benchmarks/pseudomonas_graph_baseline.json`. A relationship edge is not merge authorization. Never simplify a branch, cycle, orientation conflict, competing containment, or known forbidden edge merely to make component counts smaller.

Graph decision-policy changes must also run `pytest tests/test_decision_policy.py tests/test_decision_policy_benchmark.py` and explain every change from `benchmarks/pseudomonas_decision_policy_baseline.json`. Without explicit junction evidence, overlap components remain deferred even when their graph topology is linear. Do not silently turn graph presence into merge authorization.

Path-planning changes must run `pytest tests/test_path_planning.py tests/test_path_planning_benchmark.py`. Plans must retain every source member of every path node with explicit path-relative orientation. A path plan is metadata, not permission to construct sequence, and `benchmarks/pseudomonas_path_planning_baseline.json` must not be updated silently.

BAM/CRAM evidence changes must run `pytest tests/test_bam_evidence.py`. Tests use an injected command runner and must not require samtools. Reference names and lengths must match the sample FASTA exactly, evidence remains sample-scoped, and existing alignments must never be presented as support for a new junction.

Targeted-remapping changes must also run `pytest tests/test_junction_evidence.py`. Keep provisional-reference construction outside the evidence adapter, count distinct primary read names conservatively, retain exact samtools/minimap2 provenance, and do not connect a spanning-read count to graph eligibility without a reviewed technology-specific benchmark and explicit policy.
