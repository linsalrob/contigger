"""Contigger: conservative, provenance-aware contig reconciliation."""

from importlib.metadata import PackageNotFoundError, version

from contigger.models import RunConfig

__all__ = ["RunConfig", "__version__"]

try:
    # The distribution metadata is generated from the version in pyproject.toml.
    __version__ = version("contigger")
except PackageNotFoundError:
    # Source-tree imports before installation still need a useful, non-release value.
    __version__ = "0+unknown"
