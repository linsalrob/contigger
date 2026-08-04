"""Command-line interface for validation and conservative run planning."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from contigger import __version__
from contigger.config import build_run_config
from contigger.exceptions import ContiggerError
from contigger.manifest import ManifestValidation, parse_manifest
from contigger.merge import merge_samples
from contigger.utilities.subprocesses import find_executable


def build_parser() -> argparse.ArgumentParser:
    """Construct the public CLI parser."""
    parser = argparse.ArgumentParser(
        prog="contigger",
        description="Conservative, provenance-aware reconciliation of assembled contigs.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    validate = commands.add_parser("validate", help="validate a manifest and referenced inputs")
    validate.add_argument("--manifest", required=True, type=Path)
    validate.set_defaults(handler=_run_validate)

    merge = commands.add_parser("merge", help="plan a merge (execution is not implemented)")
    merge.add_argument("--manifest", required=True, type=Path)
    merge.add_argument("--output-prefix", required=True, type=Path)
    merge.add_argument(
        "--identity", type=float, default=98.0, help="percent identity (default: 98)"
    )
    merge.add_argument("--min-overlap", type=int, default=1000)
    merge.add_argument("--min-containment", type=int, default=500)
    merge.add_argument("--containment-coverage", type=float, default=98.0)
    merge.add_argument("--end-tolerance", type=int, default=50)
    merge.add_argument("--kmer-size", type=int, default=21)
    merge.add_argument("--window-size", type=int, default=10)
    merge.add_argument("--min-shared-minimisers", type=int, default=5)
    merge.add_argument("--max-minimiser-frequency", type=int, default=100)
    merge.add_argument("--threads", type=int, default=1)
    merge.add_argument("--evidence", choices=("none", "alignments", "reads"), default="none")
    merge.add_argument(
        "--conflict-policy",
        choices=(
            "representative",
            "majority",
            "quality-weighted",
            "sample-aware",
            "ambiguous",
            "reject",
        ),
        default="reject",
    )
    merge.add_argument("--emit-gfa", action="store_true")
    merge.add_argument("--dry-run", action="store_true")
    merge.set_defaults(handler=_run_merge)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and convert expected domain errors to exit status 2."""
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        handler = arguments.handler
        return int(handler(arguments))
    except ContiggerError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


def _run_validate(arguments: argparse.Namespace) -> int:
    validation = parse_manifest(arguments.manifest)
    _print_warnings(validation)
    print(f"validated {len(validation.samples)} sample(s)")
    return 0


def _run_merge(arguments: argparse.Namespace) -> int:
    validation = parse_manifest(arguments.manifest)
    _print_warnings(validation)
    config = build_run_config(
        identity=arguments.identity,
        min_overlap=arguments.min_overlap,
        min_containment=arguments.min_containment,
        containment_coverage=arguments.containment_coverage,
        end_tolerance=arguments.end_tolerance,
        kmer_size=arguments.kmer_size,
        window_size=arguments.window_size,
        min_shared_minimisers=arguments.min_shared_minimisers,
        max_minimiser_frequency=arguments.max_minimiser_frequency,
        threads=arguments.threads,
        evidence=arguments.evidence,
        conflict_policy=arguments.conflict_policy,
        output_prefix=arguments.output_prefix,
        emit_gfa=arguments.emit_gfa,
    )
    if arguments.dry_run:
        external_tools: dict[str, str | None] = {}
        for tool in ("minimap2", "samtools"):
            executable_path = find_executable(tool)
            external_tools[tool] = str(executable_path) if executable_path else None
        plan = {
            "configuration": config.as_dict(),
            "external_tools": external_tools,
            "mode": "dry-run",
            "samples": [
                {
                    "assembly_graph": str(sample.assembly_graph) if sample.assembly_graph else None,
                    "bam": str(sample.bam) if sample.bam else None,
                    "contigs": str(sample.contigs),
                    "sample": sample.sample,
                    "technology": sample.technology,
                }
                for sample in validation.samples
            ],
        }
        print(json.dumps(plan, indent=2, sort_keys=True))
        for tool, detected_path in external_tools.items():
            if detected_path is None:
                print(f"warning: optional external tool not found: {tool}", file=sys.stderr)
        return 0
    merge_samples(validation.samples, config)
    return 0


def _print_warnings(validation: ManifestValidation) -> None:
    for warning in validation.warnings:
        print(f"warning: {warning}", file=sys.stderr)
