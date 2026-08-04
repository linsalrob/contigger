"""Merge orchestration boundary; biological merging is deliberately absent."""

from pathlib import Path

from contigger.exceptions import FeatureNotImplementedError
from contigger.models import RunConfig, SampleInput


def merge_samples(samples: tuple[SampleInput, ...], config: RunConfig) -> tuple[Path, ...]:
    """Refuse real merging until the documented conservative pipeline exists."""
    raise FeatureNotImplementedError(
        "sequence merging is not implemented; use 'contigger merge --dry-run' to validate a plan"
    )
