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
