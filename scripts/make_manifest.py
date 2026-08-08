#!/usr/bin/env python3
"""Create a Contigger manifest from an assembly directory."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from contigger.manifest_generator import create_manifest


def main(argv: list[str] | None = None) -> int:
    """Run the standalone manifest generator."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("samples.tsv"))
    parser.add_argument("--no-recursive", action="store_true")
    arguments = parser.parse_args(argv)
    try:
        rows, warnings = create_manifest(
            arguments.directory,
            arguments.output,
            recursive=not arguments.no_recursive,
        )
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(f"wrote {len(rows)} sample(s) to {arguments.output}")
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
