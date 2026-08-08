# Contributing

Read the repository [CONTRIBUTING.md](https://github.com/linsalrob/contigger/blob/main/CONTRIBUTING.md) for the canonical contribution policy. In brief: create a focused branch, preserve conservative invariants, add deterministic tests, keep external tools behind typed adapters, and explain biological assumptions in documentation and benchmark changes.

Developer setup:

```bash
python -m pip install -e '.[dev]'
pytest
ruff check .
ruff format --check .
mypy src
```

Documentation changes should also run `python -m pip install -e '.[docs]'` and `mkdocs build --strict`.
