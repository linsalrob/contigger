# Testing

Install development dependencies with `python -m pip install -e '.[dev]'` and run:

```bash
pytest
ruff check .
ruff format --check .
mypy src
python benchmarking/validate_benchmark.py --lightweight test_data
```

External minimap2/samtools tests are marked `integration`; ordinary unit tests use injected runners or typed fake aligners. CI also runs PAF, pipeline, graph, policy, path, BAM, junction, catalogue, and candidate baselines. Documentation CI runs `python -m pip install -e '.[docs]'` followed by `mkdocs build --strict`.
