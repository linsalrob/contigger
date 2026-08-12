"""Command-line interface for validation and conservative run planning."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from contigger import __version__
from contigger.aligners.minimap2 import parse_paf
from contigger.alignment_planning import plan_selective_alignments
from contigger.benchmark import evaluate_benchmark, format_summary, write_json, write_tsv
from contigger.catalogue import (
    build_catalogue,
    catalogue_provenance,
    load_source_sequences,
    write_catalogue_fasta_path,
)
from contigger.config import build_run_config
from contigger.evidence.bam import BamEvidenceProvider
from contigger.exceptions import ContiggerError, InputValidationError
from contigger.junction_remapping_benchmark import (
    evaluate_junction_remapping_benchmark,
    format_junction_remapping_summary,
    write_junction_remapping_json,
)
from contigger.manifest import ManifestValidation, parse_manifest
from contigger.manifest_generator import create_manifest
from contigger.merge import merge_samples
from contigger.minimisers import generate_candidates, write_candidates_tsv
from contigger.models import PairRelationship
from contigger.pipeline_benchmark import (
    evaluate_pipeline_benchmark,
    format_pipeline_summary,
    write_pipeline_json,
)
from contigger.provenance import write_provenance
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.textio import open_text
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

    validate_alignments = commands.add_parser(
        "validate-alignments",
        help="validate supplied sample BAM/CRAM references against source FASTA files",
    )
    validate_alignments.add_argument("--manifest", required=True, type=Path)
    validate_alignments.add_argument("--samtools", default="samtools")
    validate_alignments.set_defaults(handler=_run_validate_alignments)

    make_manifest = commands.add_parser(
        "make-manifest", help="discover FASTA files and matching evidence sidecars"
    )
    make_manifest.add_argument("directory", type=Path)
    make_manifest.add_argument("-o", "--output", type=Path, default=Path("samples.tsv"))
    make_manifest.add_argument("--no-recursive", action="store_true")
    make_manifest.set_defaults(handler=_run_make_manifest)

    merge = commands.add_parser("merge", help="construct conservative merged contig sequences")
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
    merge.add_argument(
        "--max-candidate-pairs",
        type=int,
        help="abort before alignment if candidate pairs exceed this limit",
    )
    merge.add_argument(
        "--max-seed-pair-observations",
        type=int,
        help="abort before seed-pair expansion if its conservative upper bound exceeds this limit",
    )
    merge.add_argument(
        "--candidate-shards",
        type=int,
        default=16,
        help="disk shards for candidate evidence (default: 16; maximum: 64)",
    )
    merge.add_argument(
        "--max-queries-per-alignment-batch",
        type=int,
        default=1000,
        help="approved queries per indexed minimap2 batch (default: 1000)",
    )
    merge.add_argument("--threads", type=int, default=1)
    merge.add_argument("--minimap2-preset", choices=("asm5", "asm10", "asm20"), default="asm20")
    merge.add_argument("--index-dir", type=Path)
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

    classify_paf = commands.add_parser(
        "classify-paf", help="experimentally classify complete ordered pairs in a PAF file"
    )
    classify_paf.add_argument(
        "--paf",
        required=True,
        type=Path,
        help="input PAF (coordinates are zero-based, half-open)",
    )
    classify_paf.add_argument("--output", required=True, type=Path)
    classify_paf.add_argument("--identity", type=float, default=98.0)
    classify_paf.add_argument("--min-overlap", type=int, default=1000)
    classify_paf.add_argument("--min-containment", type=int, default=500)
    classify_paf.add_argument("--containment-coverage", type=float, default=98.0)
    classify_paf.add_argument("--end-tolerance", type=int, default=50)
    classify_paf.set_defaults(handler=_run_classify_paf)

    benchmark = commands.add_parser(
        "benchmark",
        help="score pair classification against construction-derived benchmark truth",
        description=(
            "Score ordered pairs while deferring graph-level ambiguity that requires "
            "multiple target pairs or component context. No contigs are merged."
        ),
    )
    benchmark.add_argument("--dataset", required=True, type=Path)
    benchmark.add_argument("--paf", required=True, type=Path)
    benchmark.add_argument("--output-json", type=Path)
    benchmark.add_argument("--output-tsv", type=Path)
    benchmark.add_argument("--identity", type=float, default=98.0)
    benchmark.add_argument("--min-overlap", type=int, default=1000)
    benchmark.add_argument("--min-containment", type=int, default=500)
    benchmark.add_argument("--containment-coverage", type=float, default=98.0)
    benchmark.add_argument("--end-tolerance", type=int, default=50)
    benchmark.add_argument("--fail-on-false-merge", action="store_true")
    benchmark.set_defaults(handler=_run_benchmark)

    pipeline_benchmark = commands.add_parser(
        "benchmark-pipeline",
        help="score exact catalogue, candidate recall, and eligible pair classification",
        description=(
            "Evaluate the catalogue-to-candidate-to-relationship pathway from checked-in "
            "PAF observations. This does not invoke minimap2, construct a graph, or merge contigs."
        ),
    )
    pipeline_benchmark.add_argument("--dataset", required=True, type=Path)
    pipeline_benchmark.add_argument("--paf", required=True, type=Path)
    pipeline_benchmark.add_argument("--output-json", type=Path)
    pipeline_benchmark.add_argument("--identity", type=float, default=98.0)
    pipeline_benchmark.add_argument("--min-overlap", type=int, default=1000)
    pipeline_benchmark.add_argument("--min-containment", type=int, default=500)
    pipeline_benchmark.add_argument("--containment-coverage", type=float, default=98.0)
    pipeline_benchmark.add_argument("--end-tolerance", type=int, default=50)
    pipeline_benchmark.add_argument("--kmer-size", type=int, default=21)
    pipeline_benchmark.add_argument("--window-size", type=int, default=10)
    pipeline_benchmark.add_argument("--min-shared-minimisers", type=int, default=5)
    pipeline_benchmark.add_argument("--max-minimiser-frequency", type=int, default=100)
    pipeline_benchmark.add_argument("--terminal-band", type=int, default=1000)
    pipeline_benchmark.add_argument("--fail-on-false-merge", action="store_true")
    pipeline_benchmark.set_defaults(handler=_run_pipeline_benchmark)

    junction_benchmark = commands.add_parser(
        "benchmark-junction-remapping",
        help="remap checked-in targeted FASTQs and score junction-spanning observations",
        description=(
            "Score sample-specific remapping evidence against junction truth. This command "
            "does not authorize graph edges, construct consensus, or merge contigs."
        ),
    )
    junction_benchmark.add_argument("--dataset", required=True, type=Path)
    junction_benchmark.add_argument("--minimap2", default="minimap2")
    junction_benchmark.add_argument(
        "--preset", choices=("sr", "map-ont", "map-hifi"), default="map-ont"
    )
    junction_benchmark.add_argument("--threads", type=int, default=1)
    junction_benchmark.add_argument("--minimum-spanning-flank", type=int, default=20)
    junction_benchmark.add_argument("--minimum-mapping-quality", type=int, default=0)
    junction_benchmark.add_argument("--output-json", type=Path)
    junction_benchmark.add_argument("--fail-on-false-support", action="store_true")
    junction_benchmark.set_defaults(handler=_run_junction_remapping_benchmark)

    catalogue = commands.add_parser(
        "catalogue",
        help="create a canonical exact-deduplicated sequence catalogue with provenance",
    )
    catalogue.add_argument("--manifest", required=True, type=Path)
    catalogue.add_argument("--output-fasta", required=True, type=Path)
    catalogue.add_argument("--output-provenance", required=True, type=Path)
    catalogue.set_defaults(handler=_run_catalogue)

    candidates = commands.add_parser(
        "candidates",
        help="emit positional-minimiser candidate evidence (not relationships or merges)",
    )
    candidates.add_argument("--manifest", required=True, type=Path)
    candidates.add_argument("--output", required=True, type=Path)
    candidates.add_argument("--kmer-size", type=int, default=21)
    candidates.add_argument("--window-size", type=int, default=10)
    candidates.add_argument("--min-shared-minimisers", type=int, default=5)
    candidates.add_argument("--max-minimiser-frequency", type=int, default=100)
    candidates.add_argument(
        "--terminal-band",
        type=int,
        default=1000,
        help="bases at each end eligible for terminal seed evidence (default: 1000)",
    )
    candidates.set_defaults(handler=_run_candidates)
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


def _run_validate_alignments(arguments: argparse.Namespace) -> int:
    validation = parse_manifest(arguments.manifest)
    _print_warnings(validation)
    supplied = tuple(sample for sample in validation.samples if sample.bam is not None)
    if not supplied:
        raise InputValidationError("manifest supplies no BAM/CRAM inputs")
    for sample in supplied:
        references = BamEvidenceProvider(sample, executable=arguments.samtools).validate_source()
        print(f"sample {sample.sample}: validated {len(references)} source reference(s)")
    print(f"validated BAM/CRAM inputs for {len(supplied)} sample(s)")
    return 0


def _run_make_manifest(arguments: argparse.Namespace) -> int:
    """Discover assemblies and write a valid manifest."""
    try:
        rows, warnings = create_manifest(
            arguments.directory,
            arguments.output,
            recursive=not arguments.no_recursive,
        )
    except (OSError, ValueError) as error:
        raise InputValidationError(str(error)) from error
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    print(f"wrote {len(rows)} sample(s) to {arguments.output}")
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
        minimap2_preset=arguments.minimap2_preset,
        index_dir=arguments.index_dir,
        max_candidate_pairs=arguments.max_candidate_pairs,
        max_seed_pair_observations=arguments.max_seed_pair_observations,
        candidate_shards=arguments.candidate_shards,
        max_queries_per_alignment_batch=arguments.max_queries_per_alignment_batch,
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
    outputs = merge_samples(validation.samples, config)
    print("wrote " + ", ".join(str(path) for path in outputs))
    return 0


def _print_warnings(validation: ManifestValidation) -> None:
    for warning in validation.warnings:
        print(f"warning: {warning}", file=sys.stderr)


def _run_classify_paf(arguments: argparse.Namespace) -> int:
    config = build_run_config(
        identity=arguments.identity,
        min_overlap=arguments.min_overlap,
        min_containment=arguments.min_containment,
        containment_coverage=arguments.containment_coverage,
        end_tolerance=arguments.end_tolerance,
    )
    with open_text(arguments.paf) as paf_file:
        decisions = [
            classify_pair(group, config) for group in group_ordered_pairs(parse_paf(paf_file))
        ]
    try:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("w", encoding="utf-8", newline="") as output:
            _write_relationships_tsv(decisions, output)
    except OSError as error:
        raise InputValidationError(
            f"cannot write relationships {arguments.output}: {error}"
        ) from error
    return 0


def _run_benchmark(arguments: argparse.Namespace) -> int:
    config = build_run_config(
        identity=arguments.identity,
        min_overlap=arguments.min_overlap,
        min_containment=arguments.min_containment,
        containment_coverage=arguments.containment_coverage,
        end_tolerance=arguments.end_tolerance,
    )
    report = evaluate_benchmark(arguments.dataset, arguments.paf, config)
    print(format_summary(report))
    for path, writer, label in (
        (arguments.output_json, write_json, "benchmark JSON"),
        (arguments.output_tsv, write_tsv, "benchmark TSV"),
    ):
        if path is None:
            continue
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", encoding="utf-8", newline="") as output:
                writer(report, output)
        except OSError as error:
            raise InputValidationError(f"cannot write {label} {path}: {error}") from error
    return int(arguments.fail_on_false_merge and report.summary.false_merges > 0)


def _run_pipeline_benchmark(arguments: argparse.Namespace) -> int:
    config = build_run_config(
        identity=arguments.identity,
        min_overlap=arguments.min_overlap,
        min_containment=arguments.min_containment,
        containment_coverage=arguments.containment_coverage,
        end_tolerance=arguments.end_tolerance,
    )
    report = evaluate_pipeline_benchmark(
        arguments.dataset,
        arguments.paf,
        config,
        kmer_size=arguments.kmer_size,
        window_size=arguments.window_size,
        min_shared_minimisers=arguments.min_shared_minimisers,
        max_minimiser_frequency=arguments.max_minimiser_frequency,
        terminal_band=arguments.terminal_band,
    )
    print(format_pipeline_summary(report))
    if arguments.output_json is not None:
        try:
            arguments.output_json.parent.mkdir(parents=True, exist_ok=True)
            with arguments.output_json.open("w", encoding="utf-8", newline="") as output:
                write_pipeline_json(report, output)
        except OSError as error:
            raise InputValidationError(
                f"cannot write pipeline benchmark JSON {arguments.output_json}: {error}"
            ) from error
    return int(arguments.fail_on_false_merge and report.summary.relationship_stage_false_merges > 0)


def _run_junction_remapping_benchmark(arguments: argparse.Namespace) -> int:
    report = evaluate_junction_remapping_benchmark(
        arguments.dataset,
        minimap2=arguments.minimap2,
        preset=arguments.preset,
        threads=arguments.threads,
        minimum_spanning_flank=arguments.minimum_spanning_flank,
        minimum_mapping_quality=arguments.minimum_mapping_quality,
    )
    print(format_junction_remapping_summary(report))
    if arguments.output_json is not None:
        try:
            arguments.output_json.parent.mkdir(parents=True, exist_ok=True)
            with arguments.output_json.open("w", encoding="utf-8", newline="") as output:
                write_junction_remapping_json(report, output)
        except OSError as error:
            raise InputValidationError(
                f"cannot write junction-remapping benchmark JSON {arguments.output_json}: {error}"
            ) from error
    return int(arguments.fail_on_false_support and report.summary.false_support_sample_cases > 0)


def _run_catalogue(arguments: argparse.Namespace) -> int:
    validation = parse_manifest(arguments.manifest)
    _print_warnings(validation)
    catalogue = build_catalogue(load_source_sequences(validation.samples))
    write_catalogue_fasta_path(catalogue, arguments.output_fasta)
    try:
        arguments.output_provenance.parent.mkdir(parents=True, exist_ok=True)
        write_provenance(arguments.output_provenance, catalogue_provenance(catalogue))
    except OSError as error:
        raise InputValidationError(
            f"cannot write catalogue provenance {arguments.output_provenance}: {error}"
        ) from error
    print(
        f"catalogued {len(catalogue.members)} source contig(s) as "
        f"{len(catalogue.sequences)} canonical sequence(s)"
    )
    return 0


def _run_candidates(arguments: argparse.Namespace) -> int:
    validation = parse_manifest(arguments.manifest)
    _print_warnings(validation)
    catalogue = build_catalogue(load_source_sequences(validation.samples))
    candidates = generate_candidates(
        catalogue.sequences,
        kmer_size=arguments.kmer_size,
        window_size=arguments.window_size,
        min_shared_minimisers=arguments.min_shared_minimisers,
        max_minimiser_frequency=arguments.max_minimiser_frequency,
        terminal_band=arguments.terminal_band,
    )
    requests = plan_selective_alignments(catalogue.sequences, candidates)
    try:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("w", encoding="utf-8", newline="") as output:
            write_candidates_tsv(candidates, output)
    except OSError as error:
        raise InputValidationError(
            f"cannot write candidates {arguments.output}: {error}"
        ) from error
    print(
        f"planned {len(requests)} selective alignment candidate(s) from "
        f"{len(catalogue.sequences)} canonical sequence(s)"
    )
    return 0


def _write_relationships_tsv(decisions: list[PairRelationship], output: TextIO) -> None:
    """Write deterministic diagnostic relationships with zero-based half-open coordinates."""
    columns = (
        "query",
        "target",
        "relationship_type",
        "orientation",
        "identity",
        "aligned_length",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "query_coverage",
        "target_coverage",
        "status",
        "accepted_hits",
        "rejected_hits",
        "reasons",
    )
    output.write("\t".join(columns) + "\n")
    for decision in decisions:
        hit = decision.representative_hit
        if hit is None and decision.accepted_hits:
            hit = decision.accepted_hits[0]
        if hit is None:
            hit = decision.rejected_alignments[0].hit
        relationship = decision.relationship
        reasons = tuple(
            sorted(
                set(relationship.reasons)
                | {
                    reason
                    for rejected in decision.rejected_alignments
                    for reason in rejected.relationship.reasons
                }
            )
        )
        row = (
            hit.query_id,
            hit.target_id,
            relationship.relationship_type.value,
            relationship.orientation.value,
            f"{relationship.identity:.6f}",
            str(relationship.aligned_length),
            str(hit.query_start),
            str(hit.query_end),
            str(hit.target_start),
            str(hit.target_end),
            f"{relationship.query_coverage:.6f}",
            f"{relationship.target_coverage:.6f}",
            relationship.status,
            str(len(decision.accepted_hits)),
            str(len(decision.rejected_alignments)),
            "; ".join(reasons),
        )
        output.write("\t".join(row) + "\n")
