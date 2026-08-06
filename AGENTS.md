# Instructions for coding agents

## Project summary

Contigger conservatively reconciles assembled contigs across samples while retaining complete provenance. It is intended for metagenomic, microbial, viral, and phage contigs, where repeats, strain variation, assembly errors, and uneven evidence make aggressive joins dangerous.

## Primary design objective

> A missed merge is preferable to a false merge.

## Non-negotiable invariants

1. Never merge contigs solely because they exceed a sequence identity threshold.
2. A mergeable overlap must have compatible terminal geometry.
3. Internal local similarity is not sufficient evidence for a merge.
4. Containment and terminal overlap are different relationship types.
5. Reverse-complement orientation must always be explicit.
6. Ambiguous graph branches must be preserved rather than greedily collapsed.
7. Sample-specific read evidence must not be blindly pooled.
8. Existing BAM alignments can support bases within source contigs but cannot directly prove a newly constructed junction.
9. New junction validation may require targeted remapping of relevant reads.
10. Every removed, contained, or merged contig must remain recoverable through provenance.
11. Output ordering must be deterministic.
12. Coordinates must use one documented convention internally.
13. No command may silently ignore malformed input, missing references, identifier collisions, or inconsistent sequence lengths.
14. No placeholder implementation may report successful biological results.

## Coordinate convention

Internal coordinates are zero-based and half-open. External formats must be converted explicitly at input and output boundaries. Strand or orientation must never be inferred from coordinate order alone; use an explicit `Orientation` value.

## Development workflow

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
pytest
ruff check .
ruff format --check .
mypy src
```

Combined verification:

```bash
pytest && ruff check . && ruff format --check . && mypy src
```

## Coding standards

- Give every function a typed signature and every public class or function a docstring.
- Keep modules and functions small and focused; use dataclasses or similarly explicit typed models.
- Raise clear domain exceptions rather than generic `ValueError` where domain meaning matters.
- Never use `shell=True`. Pass external commands as argument arrays, record exact commands, capture tool versions, and include stderr on failure.
- Do not hide configuration in global state. Validate explicitly at every public interface.
- Sort deterministically before writing output; never rely on dictionary or graph traversal order.
- Preserve public data models so a native or Rust backend can replace Python internals without unnecessary CLI changes.
- Comments should explain biological or algorithmic reasoning, not restate syntax.

## Testing standards

- Unit tests must not require minimap2, samtools, or large datasets.
- Mark integration tests requiring external tools with `pytest.mark.integration`.
- Synthetic fixtures must cover orientation, containment, overlap geometry, ambiguity, and identity thresholds.
- Test classification first with manually constructed alignment records, not an invoked aligner.
- Every regression involving a false join receives a permanent test case.
- Changes affecting relationship classification must run both checked-in Pseudomonas PAF benchmarks.
- Benchmark baselines must never be updated silently; explain every count change and add a permanent regression test for each new false merge.
- Catalogue or candidate changes must run the exact/RC and Pseudomonas candidate baselines; candidate recall losses require explicit review and must not be hidden by baseline updates.
- Catalogue, candidate, alignment-planning, or relationship changes must run both checked-in pipeline benchmarks; staged baseline changes require an explanation and must never be updated silently.
- Graph model or construction changes must run the checked-in graph regressions; never hide a forbidden edge or ambiguity loss by silently updating the graph baseline.
- Graph decision-policy changes must run the checked-in policy regression; no overlap may become eligible without explicit junction evidence, and policy baselines must never hide a forbidden edge.
- Path-planning changes must retain every catalogue source member with explicit orientation, preserve deferred components, and run the checked-in path-planning regression without emitting sequence.
- Use fixed seeds in randomised tests and pytest temporary directories for temporary files.

## External tool policy

The initial implementation may wrap minimap2 for sequence alignment and targeted read remapping, and samtools for BAM/CRAM validation, indexing, extraction, and conversion. Each wrapper must check executable availability, capture the exact command and tool version, report stderr on failure, avoid shell interpolation, and be replaceable by a mock implementation in tests.

## Definition of done

A feature is complete only when its implementation is typed; tests cover expected and failure behaviour; documentation is updated; output remains deterministic; provenance requirements are satisfied; biological assumptions are explicit; and tests, linting, formatting, and type checking pass.

## Agent behaviour

- Inspect this file, `DESIGN.md`, and relevant tests before changing code.
- Make narrow changes and avoid rewriting unrelated modules.
- Preserve public interfaces unless a change is explicitly justified.
- Add tests before or alongside algorithmic changes.
- State uncertainties instead of inventing biological rules.
- Leave explicit TODOs when a decision requires benchmarking or domain review.
- Never replace conservative logic with a more aggressive heuristic without documentation and regression tests.
