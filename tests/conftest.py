"""Shared pytest configuration for source-layout tests."""

import sys
from pathlib import Path

SRC = Path(__file__).parents[1] / "src"
sys.path.insert(0, str(SRC))
