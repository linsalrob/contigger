#!/usr/bin/env python3
"""Build the deterministic Contigger validation benchmark (truth is construction-derived)."""

import argparse
import collections
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from benchmarklib.fasta import format_fasta, read_fasta, reverse_complement  # noqa: E402
from benchmarklib.intervals import extract  # noqa: E402
from benchmarklib.mutations import delete, insert, substitute  # noqa: E402
from benchmarklib.reads import partition, stable_hash  # noqa: E402
from benchmarklib.truth import sort_relationships  # noqa: E402
from benchmarklib.utilities import (  # noqa: E402
    printable,
    sha256,
    write_gzip,
    write_json,
    write_tsv,
)


def tool(name):
    configured = os.environ.get(name.upper())
    resolved = shutil.which(configured or name)
    if resolved:
        return resolved
    detail = f" configured by {name.upper()}={configured!r}" if configured else " on PATH"
    raise SystemExit(f"ERROR: required installed tool not found{detail}: {name}")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=str(HERE / "config.yaml"))
    p.add_argument("--output", default="test_data")
    p.add_argument("--reads")
    p.add_argument("--assembly")
    p.add_argument("--graph")
    p.add_argument("--seed", type=int)
    p.add_argument("--force", action="store_true")
    p.add_argument("--keep-work", action="store_true")
    return p.parse_args()


def source_stats(reads, assembly, graph, seed):
    for p in (reads, assembly, graph):
        if not p.is_file() or not os.access(p, os.R_OK):
            raise SystemExit(f"ERROR: missing/unreadable input: {p}")
    recs = read_fasta(assembly)
    if len(recs) != 1 or not recs[0][1]:
        raise SystemExit("ERROR: assembly must contain exactly one non-empty sequence")
    aid, seq = recs[0]
    header = assembly.open().readline().strip()
    n = len(seq)
    if not 6_000_000 < n < 8_000_000:
        raise SystemExit(f"ERROR: materially inconsistent assembly length: {n}")
    g = collections.Counter()
    seg_lengths = []
    with graph.open() as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip().split("\t")
                g["segments"] += 1
                seg_lengths.append(len(f[2]) if f[2] != "*" else 0)
            elif line.startswith("L\t"):
                f = line.rstrip().split("\t")
                g["links"] += 1
                g["self_links"] += f[1] == f[3]
            elif line.startswith("P\t"):
                g["paths"] += 1
            elif line.startswith("W\t"):
                g["walks"] += 1
    gfa_ok = g["segments"] == 1 and seg_lengths == [n] and g["self_links"] >= 1
    if not gfa_ok:
        raise SystemExit("ERROR: GFA materially inconsistent with one circular assembly")
    # Source statistics are independently streamed with seqkit during build verification.
    cmd = [tool("seqkit"), "stats", "-a", "-T", str(reads)]
    text = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL)
    line = text.strip().splitlines()[-1].split("\t")
    hdr = text.strip().splitlines()[-2].split("\t")
    st = dict(zip(hdr, line, strict=True))
    rc = int(st["num_seqs"])
    bases = int(st["sum_len"])
    if abs(rc - 2_009_823) > 20_000 or abs(bases - 6_229_995_885) > 100_000_000:
        raise SystemExit("ERROR: reads materially inconsistent with expected source")
    return seq, {
        "source_basenames": {"reads": reads.name, "assembly": assembly.name, "graph": graph.name},
        "source_sha256": {
            "reads": sha256(reads),
            "assembly": sha256(assembly),
            "graph": sha256(graph),
        },
        "source_file_sizes": {
            "reads": reads.stat().st_size,
            "assembly": assembly.stat().st_size,
            "graph": graph.stat().st_size,
        },
        "reads": {
            "count": rc,
            "total_bases": bases,
            "mean_length": float(st["avg_len"]),
            "median_length": int(st["Q2"]),
            "n50": int(st["N50"]),
            "minimum_length": int(st["min_len"]),
            "maximum_length": int(st["max_len"]),
            "gc_percent": float(st["GC(%)"]),
        },
        "assembly": {
            "contig_id": aid,
            "length": n,
            "gc_percent": 100 * (seq.upper().count("G") + seq.upper().count("C")) / n,
            "header": header,
            "header_claims_circular": (
                "circular=true" in header.lower() or "topology=circular" in header.lower()
            ),
        },
        "gfa": {
            "segments": g["segments"],
            "segment_lengths": seg_lengths,
            "links": g["links"],
            "self_links": g["self_links"],
            "paths": g["paths"],
            "walks": g["walks"],
            "consistent_single_circular": gfa_ok,
        },
        "inferred_circularity": gfa_ok and ("circular" in header.lower()),
        "generation_date": "2026-08-05",
        "benchmark_seed": seed,
    }


def build_sequences(seq, cfg):
    seed = cfg["seed"]
    contigs = {}
    meta = []
    rel = []
    merges = {}
    muts = []
    intervals = []
    cases = []
    junctions = []
    segments = []
    samples = ["S01", "S02", "S03"]

    def add(
        cid,
        sample,
        case,
        s,
        e,
        st="+",
        wrap=False,
        sequence=None,
        ctype="reference_slice",
        desc="",
        mutset="",
        segs=None,
    ):
        x = sequence if sequence is not None else extract(seq, s, e, st, wrap)
        contigs[cid] = (sample, x)
        meta.append(
            dict(
                contig_id=cid,
                sample=sample,
                case_id=case,
                source_start=(s if segs is None else ""),
                source_end=(e if segs is None else ""),
                strand=(st if segs is None else ""),
                circular_wrap=str(wrap).lower(),
                length=len(x),
                mutation_set=mutset,
                construction_type=ctype,
                description=desc,
            )
        )
        if segs:
            for i, (a, b, z) in enumerate(segs, 1):
                segments.append(
                    dict(contig_id=cid, segment_index=i, source_start=a, source_end=b, strand=z)
                )
        intervals.append(
            dict(
                case_id=case,
                contig_id=cid,
                source_start=(s if segs is None else ""),
                source_end=(e if segs is None else ""),
                strand=st,
                circular_wrap=str(wrap).lower(),
                selection_rationale=desc,
            )
        )
        return x

    def relationship(
        case, q, t, kind, orient, status, allowed, ov, ident, qs, qe, ts, te, reason, amb=""
    ):
        rel.append(
            dict(
                query_id=q,
                target_id=t,
                expected_relationship=kind,
                expected_orientation=orient,
                expected_status=status,
                merge_allowed=str(allowed).lower(),
                expected_overlap_length=ov,
                expected_identity=f"{ident:.4f}",
                query_start=qs,
                query_end=qe,
                target_start=ts,
                target_end=te,
                ambiguity_group=amb,
                case_id=case,
                reason=reason,
            )
        )

    def case(cid, category, desc):
        cases.append(
            dict(
                case_id=cid,
                category=category,
                description=desc,
                governing_decision="merge only when uniquely and conservatively supported",
            )
        )

    def overlap_pair(
        caseid,
        leftid,
        rightid,
        start,
        ov=5000,
        leftlen=20000,
        rightlen=20000,
        mut=None,
        category="positive",
        reason="constructed native suffix-prefix overlap",
        end_extra=0,
    ):
        case(caseid, category, reason)
        junction = start + leftlen
        L = add(
            leftid,
            samples[len(contigs) % 3],
            caseid,
            start,
            start + leftlen + end_extra,
            desc="unique-locus construction; deterministic spaced locus",
        )
        R = extract(seq, junction - ov, junction - ov + rightlen)
        mutset = ""
        target_overlap_span = ov
        if mut:
            typ, arg = mut
            if typ == "sub":
                portion, newrows = R[:ov], None
                changed, newrows = substitute(portion, arg, seed + arg)
                R = changed + R[ov:]
                mutset = f"{caseid}_substitutions"
                for p, a, b in newrows:
                    muts.append(
                        dict(
                            mutation_set=mutset,
                            contig_id=rightid,
                            position=p,
                            mutation_type="substitution",
                            reference=a,
                            alternate=b,
                            case_id=caseid,
                        )
                    )
            elif typ == "ins":
                R = insert(R, ov // 2, arg)
                target_overlap_span += len(arg)
                mutset = f"{caseid}_insertion"
                muts.append(
                    dict(
                        mutation_set=mutset,
                        contig_id=rightid,
                        position=ov // 2,
                        mutation_type="insertion",
                        reference="-",
                        alternate=arg,
                        case_id=caseid,
                    )
                )
            elif typ == "del":
                R = delete(R, ov // 2, arg)
                target_overlap_span -= arg
                mutset = f"{caseid}_deletion"
                muts.append(
                    dict(
                        mutation_set=mutset,
                        contig_id=rightid,
                        position=ov // 2,
                        mutation_type="deletion",
                        reference=extract(
                            seq, junction - ov + ov // 2, junction - ov + ov // 2 + arg
                        ),
                        alternate="-",
                        case_id=caseid,
                    )
                )
        add(
            rightid,
            samples[len(contigs) % 3],
            caseid,
            junction - ov,
            junction - ov + rightlen,
            sequence=R,
            ctype="mutated_slice" if mut else "reference_slice",
            desc=reason,
            mutset=mutset,
        )
        ident = 100.0 if not mut else (100 - 100 * (arg if isinstance(arg, int) else len(arg)) / ov)
        allowed = (
            category == "positive"
            and ov >= cfg["thresholds"]["minimum_overlap"]
            and ident >= cfg["thresholds"]["identity_percent"]
            and end_extra <= cfg["thresholds"]["end_tolerance"]
        )
        status = "VALID" if allowed else "FORBIDDEN"
        relationship(
            caseid,
            leftid,
            rightid,
            "QUERY_SUFFIX_TO_TARGET_PREFIX",
            "+",
            status,
            allowed,
            ov,
            ident,
            len(L) - ov - end_extra,
            len(L) - end_extra,
            0,
            target_overlap_span,
            reason,
        )
        relationship(
            caseid,
            rightid,
            leftid,
            "TARGET_SUFFIX_TO_QUERY_PREFIX",
            "+",
            status,
            allowed,
            ov,
            ident,
            0,
            target_overlap_span,
            len(L) - ov - end_extra,
            len(L) - end_extra,
            reason,
        )
        if allowed:
            # The left contig supplies the overlap. For indels, consume the target's
            # actual aligned span before appending its unique suffix.
            merges[caseid] = L + R[target_overlap_span:]
        junctions.append(
            dict(
                case_id=caseid,
                left_contig=leftid,
                right_contig=rightid,
                junction_is_true=str(category == "positive").lower(),
                source_coordinate=junction,
                circular_wrap="false",
                selected_spanning_reads=0,
                selected_nonspanning_reads=0,
                expected_read_evidence="native spanning reads requested"
                if category == "positive"
                else "no read should support artificial adjacency",
                reason=reason,
            )
        )

    pos = 100000
    for r in range(1, 4):
        tag = f"r{r:02d}"
        # exact and reverse exact
        for rev in (False, True):
            cid = f"{'reverse_exact' if rev else 'exact'}_{tag}"
            case(
                cid,
                "positive",
                "exact sequence equality" + (" under reverse complement" if rev else ""),
            )
            base = extract(seq, pos, pos + 20000)
            add(
                cid + "_a",
                samples[(r - 1) % 3],
                cid,
                pos,
                pos + 20000,
                desc="unique-locus exact duplicate",
            )
            add(
                cid + "_b",
                samples[r % 3],
                cid,
                pos,
                pos + 20000,
                "-" if rev else "+",
                sequence=reverse_complement(base) if rev else base,
                desc="reverse-complement duplicate" if rev else "exact duplicate",
            )
            relationship(
                cid,
                cid + "_a",
                cid + "_b",
                "EXACT_MATCH",
                "-" if rev else "+",
                "VALID",
                True,
                20000,
                100,
                0,
                20000,
                0,
                20000,
                "exact construction",
            )
            relationship(
                cid,
                cid + "_b",
                cid + "_a",
                "EXACT_MATCH",
                "-" if rev else "+",
                "VALID",
                True,
                20000,
                100,
                0,
                20000,
                0,
                20000,
                "exact construction",
            )
            pos += 30000
        # both containment directions represented by ordered rows
        cc = f"contain_{tag}"
        case(cc, "positive", "10 kb query exactly contained in 25 kb target")
        add(
            f"contain_{tag}_long",
            samples[r % 3],
            cc,
            pos,
            pos + 25000,
            desc="unique containment target",
        )
        add(
            f"contain_{tag}_short",
            samples[(r + 1) % 3],
            cc,
            pos + 7000,
            pos + 17000,
            desc="exact contained slice",
        )
        relationship(
            cc,
            f"contain_{tag}_short",
            f"contain_{tag}_long",
            "QUERY_CONTAINED_IN_TARGET",
            "+",
            "VALID",
            True,
            10000,
            100,
            0,
            10000,
            7000,
            17000,
            "constructed containment",
        )
        relationship(
            cc,
            f"contain_{tag}_long",
            f"contain_{tag}_short",
            "TARGET_CONTAINED_IN_QUERY",
            "+",
            "VALID",
            True,
            10000,
            100,
            7000,
            17000,
            0,
            10000,
            "ordered inverse containment",
        )
        pos += 40000
        overlap_pair(
            f"overlap_forward_{tag}",
            f"overlap_forward_{tag}_left",
            f"overlap_forward_{tag}_right",
            pos,
        )
        pos += 50000
        overlap_pair(
            f"overlap_reverse_direction_{tag}",
            f"overlap_reverse_direction_{tag}_left",
            f"overlap_reverse_direction_{tag}_right",
            pos,
        )
        pos += 50000
        # RC terminal overlap
        cid = f"overlap_rc_{tag}"
        case(cid, "positive", "terminal overlap requires reverse-complement target")
        L = add(
            cid + "_left", samples[r % 3], cid, pos, pos + 20000, desc="unique RC overlap locus"
        )
        native = extract(seq, pos + 15000, pos + 35000)
        R = reverse_complement(native)
        add(
            cid + "_right",
            samples[(r + 1) % 3],
            cid,
            pos + 15000,
            pos + 35000,
            "-",
            sequence=R,
            desc="reverse-complement terminal partner",
        )
        relationship(
            cid,
            cid + "_left",
            cid + "_right",
            "QUERY_SUFFIX_TO_TARGET_PREFIX",
            "-",
            "VALID",
            True,
            5000,
            100,
            15000,
            20000,
            15000,
            20000,
            "reverse complement restores suffix-prefix overlap",
        )
        merges[cid] = L + native[5000:]
        pos += 50000
    overlap_pair(
        "small_insertion",
        "small_insertion_left",
        "small_insertion_right",
        pos,
        mut=("ins", "ACG"),
        reason="3 bp insertion inside terminal overlap; source reads support reference allele",
    )
    pos += 50000
    overlap_pair(
        "small_deletion",
        "small_deletion_left",
        "small_deletion_right",
        pos,
        mut=("del", 3),
        reason="3 bp deletion inside terminal overlap; source reads support reference allele",
    )
    pos += 50000
    for pct, nsub in (("9801", 199), ("9800", 200), ("9799", 201)):
        overlap_pair(
            f"identity_{pct}",
            f"identity_{pct}_left",
            f"identity_{pct}_right",
            pos,
            ov=10000,
            mut=("sub", nsub),
            category="positive" if nsub <= 200 else "negative",
            reason=f"10 kb substitution-only boundary: {nsub} substitutions",
        )
        pos += 50000
    for ov in (999, 1000, 1001):
        overlap_pair(
            f"length_{ov}",
            f"length_{ov}_left",
            f"length_{ov}_right",
            pos,
            ov=ov,
            category="positive" if ov >= 1000 else "negative",
            reason=f"minimum-overlap boundary at {ov} bp",
        )
        pos += 50000
    for gap in (49, 50, 51):
        overlap_pair(
            f"end_tolerance_{gap}",
            f"end_tolerance_{gap}_left",
            f"end_tolerance_{gap}_right",
            pos,
            end_extra=gap,
            category="positive" if gap <= 50 else "negative",
            reason=f"alignment endpoint is {gap} bp from query end",
        )
        pos += 50000
    # containment coverage boundary: 10 kb query with exact aligned fractions
    for label, aligned in (("below", 9799), ("at", 9800), ("above", 9801)):
        cid = f"containment_coverage_{label}"
        case(
            cid,
            "positive" if aligned >= 9800 else "negative",
            f"containment aligned coverage {aligned / 100:.2f}%",
        )
        add(
            cid + "_long",
            samples[len(contigs) % 3],
            cid,
            pos,
            pos + 25000,
            desc="containment boundary target",
        )
        short_seq = extract(seq, pos + 5000, pos + 5000 + aligned) + extract(
            seq, pos + 40000, pos + 40000 + (10000 - aligned)
        )
        add(
            cid + "_short",
            samples[len(contigs) % 3],
            cid,
            0,
            0,
            sequence=short_seq,
            ctype="chimeric_tail",
            desc=f"{aligned} aligned and {10000 - aligned} unaligned bases",
            segs=[
                (pos + 5000, pos + 5000 + aligned, "+"),
                (pos + 40000, pos + 40000 + 10000 - aligned, "+"),
            ],
        )
        allowed = aligned >= 9800
        relationship(
            cid,
            cid + "_short",
            cid + "_long",
            "QUERY_CONTAINED_IN_TARGET" if allowed else "NO_RELATIONSHIP",
            "+",
            "VALID" if allowed else "FORBIDDEN",
            allowed,
            aligned,
            100,
            0,
            aligned,
            5000,
            5000 + aligned,
            f"coverage={aligned / 100:.2f}%; inclusive 98% threshold",
        )
        pos += 60000
    # Circular origin
    cid = "circular_origin"
    case(cid, "positive", "native overlap crosses reference end/origin")
    n = len(seq)
    L = add(
        cid + "_left",
        "S01",
        cid,
        n - 15000,
        n + 5000,
        "+",
        True,
        desc="circular slice crossing origin",
    )
    R = add(cid + "_right", "S03", cid, 0, 20000, desc="origin-proximal slice")
    relationship(
        cid,
        cid + "_left",
        cid + "_right",
        "QUERY_SUFFIX_TO_TARGET_PREFIX",
        "+",
        "VALID",
        True,
        5000,
        100,
        15000,
        20000,
        0,
        5000,
        "circular native junction",
    )
    merges[cid] = L + R[5000:]
    junctions.append(
        dict(
            case_id=cid,
            left_contig=cid + "_left",
            right_contig=cid + "_right",
            junction_is_true="true",
            source_coordinate=0,
            circular_wrap="true",
            selected_spanning_reads=0,
            selected_nonspanning_reads=0,
            expected_read_evidence="reads spanning origin where available",
            reason="native circular adjacency",
        )
    )
    # negatives, deliberately conservative
    negatives = [
        ("unrelated_unique", "two distant unique contigs have no constructed relationship"),
        ("internal_shared", "shared high-identity block is internal to both contigs"),
        ("repeat_ambiguity", "synthetic repeat has two incompatible terminal placements"),
        (
            "incompatible_placements",
            "two valid pairwise placements cannot both define one adjacency",
        ),
        ("opposite_orientations", "non-equivalent alignments support opposite orientations"),
        ("chimera", "contig joins two distant source loci"),
        ("low_complexity", "homopolymer-rich terminal match is non-specific"),
        ("terminal_repeat", "terminal synthetic repeat lacks unique-flank adjacency evidence"),
        ("near_containment_fail", "near containment fails end/coverage criterion"),
        ("high_identity_short", "high identity overlap is only 900 bp"),
    ]
    repeat = extract(seq, 2_000_000, 2_006_000)
    for i, (cid, reason) in enumerate(negatives):
        case(
            cid,
            "ambiguous"
            if cid
            in {
                "repeat_ambiguity",
                "incompatible_placements",
                "opposite_orientations",
                "terminal_repeat",
            }
            else "negative",
            reason,
        )
        a0 = 2_200_000 + i * 80000
        b0 = 4_000_000 + i * 70000
        if cid in {
            "repeat_ambiguity",
            "terminal_repeat",
            "incompatible_placements",
            "opposite_orientations",
        }:
            A = extract(seq, a0, a0 + 14000) + repeat
            B = repeat + extract(seq, b0, b0 + 14000)
            C = repeat + extract(seq, b0 + 30000, b0 + 44000)
            add(
                cid + "_left",
                "S01",
                cid,
                0,
                0,
                sequence=A,
                ctype="synthetic_repeat",
                desc=reason,
                segs=[(a0, a0 + 14000, "+"), (2_000_000, 2_006_000, "+")],
            )
            add(
                cid + "_right_a",
                "S02",
                cid,
                0,
                0,
                sequence=B,
                ctype="synthetic_repeat",
                desc=reason,
                segs=[(2_000_000, 2_006_000, "+"), (b0, b0 + 14000, "+")],
            )
            add(
                cid + "_right_b",
                "S03",
                cid,
                0,
                0,
                sequence=C,
                ctype="synthetic_repeat",
                desc=reason,
                segs=[(2_000_000, 2_006_000, "+"), (b0 + 30000, b0 + 44000, "+")],
            )
            for t in (cid + "_right_a", cid + "_right_b"):
                relationship(
                    cid,
                    cid + "_left",
                    t,
                    "AMBIGUOUS_OVERLAP",
                    "+",
                    "AMBIGUOUS",
                    False,
                    6000,
                    100,
                    14000,
                    20000,
                    0,
                    6000,
                    reason,
                    cid,
                )
        elif cid == "internal_shared":
            shared = extract(seq, a0, a0 + 5000)
            A = extract(seq, a0 - 7000, a0) + shared + extract(seq, a0 + 5000, a0 + 12000)
            B = extract(seq, b0, b0 + 7000) + shared + extract(seq, b0 + 7000, b0 + 14000)
            add(
                cid + "_a",
                "S01",
                cid,
                0,
                0,
                sequence=A,
                ctype="synthetic_internal_match",
                desc=reason,
                segs=[(a0 - 7000, a0 + 12000, "+")],
            )
            add(
                cid + "_b",
                "S02",
                cid,
                0,
                0,
                sequence=B,
                ctype="synthetic_internal_match",
                desc=reason,
                segs=[(b0, b0 + 7000, "+"), (a0, a0 + 5000, "+"), (b0 + 7000, b0 + 14000, "+")],
            )
            relationship(
                cid,
                cid + "_a",
                cid + "_b",
                "NO_RELATIONSHIP",
                "+",
                "FORBIDDEN",
                False,
                5000,
                100,
                7000,
                12000,
                7000,
                12000,
                reason,
            )
        elif cid == "chimera":
            A = extract(seq, a0, a0 + 10000) + extract(seq, b0, b0 + 10000)
            add(
                cid + "_contig",
                "S01",
                cid,
                0,
                0,
                sequence=A,
                ctype="chimera",
                desc=reason,
                segs=[(a0, a0 + 10000, "+"), (b0, b0 + 10000, "+")],
            )
            B = add(cid + "_native", "S03", cid, a0 - 5000, a0 + 15000, desc="native comparator")
            relationship(
                cid,
                cid + "_contig",
                cid + "_native",
                "NO_RELATIONSHIP",
                "+",
                "FORBIDDEN",
                False,
                10000,
                100,
                0,
                10000,
                5000,
                15000,
                reason,
            )
        elif cid == "low_complexity":
            A = extract(seq, a0, a0 + 15000) + "A" * 5000
            B = "A" * 5000 + extract(seq, b0, b0 + 15000)
            add(
                cid + "_a",
                "S01",
                cid,
                0,
                0,
                sequence=A,
                ctype="synthetic_low_complexity",
                desc=reason,
                segs=[(a0, a0 + 15000, "+")],
            )
            add(
                cid + "_b",
                "S02",
                cid,
                0,
                0,
                sequence=B,
                ctype="synthetic_low_complexity",
                desc=reason,
                segs=[(b0, b0 + 15000, "+")],
            )
            relationship(
                cid,
                cid + "_a",
                cid + "_b",
                "NO_RELATIONSHIP",
                "+",
                "FORBIDDEN",
                False,
                5000,
                100,
                15000,
                20000,
                0,
                5000,
                reason,
            )
        elif cid == "high_identity_short":
            overlap_pair(
                cid, cid + "_left", cid + "_right", a0, ov=900, category="negative", reason=reason
            )
        elif cid == "near_containment_fail":
            L = add(cid + "_long", "S02", cid, a0, a0 + 25000, desc=reason)
            S = extract(seq, a0 + 500, a0 + 10300) + extract(seq, b0, b0 + 200)
            add(
                cid + "_short",
                "S03",
                cid,
                0,
                0,
                sequence=S,
                ctype="chimeric_tail",
                desc=reason,
                segs=[(a0 + 500, a0 + 10300, "+"), (b0, b0 + 200, "+")],
            )
            relationship(
                cid,
                cid + "_short",
                cid + "_long",
                "NO_RELATIONSHIP",
                "+",
                "FORBIDDEN",
                False,
                9800,
                100,
                0,
                9800,
                500,
                10300,
                reason,
            )
        else:
            add(cid + "_a", "S01", cid, a0, a0 + 20000, desc=reason)
            add(cid + "_b", "S03", cid, b0, b0 + 20000, desc=reason)
            relationship(
                cid,
                cid + "_a",
                cid + "_b",
                "NO_RELATIONSHIP",
                ".",
                "FORBIDDEN",
                False,
                0,
                0,
                0,
                0,
                0,
                0,
                reason,
            )
    return (
        contigs,
        meta,
        sort_relationships(rel),
        merges,
        muts,
        intervals,
        cases,
        junctions,
        segments,
    )


def main():
    a = parse_args()
    cfg = yaml.safe_load(Path(a.config).read_text())
    cfg["seed"] = a.seed if a.seed is not None else cfg["seed"]
    for k, v in (("reads", a.reads), ("assembly", a.assembly), ("graph", a.graph)):
        if v:
            cfg["inputs"][k] = v
    paths = {
        k: (ROOT / v).resolve() if not Path(v).is_absolute() else Path(v)
        for k, v in cfg["inputs"].items()
    }
    out = Path(a.output).resolve()
    work = HERE / ".work"
    if out.exists():
        if not a.force:
            raise SystemExit(f"ERROR: output exists: {out}; use --force")
        shutil.rmtree(out)
    out.mkdir()
    [(out / d).mkdir() for d in ("metadata", "contigs", "reads", "alignments", "expected", "cases")]
    work.mkdir(exist_ok=True)
    seq, source = source_stats(paths["reads"], paths["assembly"], paths["graph"], cfg["seed"])
    commands = []
    # Required private index/self-alignment/whole-genome BAM; reuse valid intermediates.
    mmi = work / "source.mmi"
    wbam = work / "source_reads.bam"
    if not mmi.exists():
        subprocess.run([tool("minimap2"), "-d", mmi, paths["assembly"]], check=True)
        commands.append(printable(["minimap2", "-d", mmi.name, paths["assembly"].name]))
    selfp = work / "source_self.asm5.paf"
    if not selfp.exists():
        with selfp.open("wb") as f:
            subprocess.run(
                [
                    tool("minimap2"),
                    "-cx",
                    "asm5",
                    "--secondary=yes",
                    paths["assembly"],
                    paths["assembly"],
                ],
                check=True,
                stdout=f,
            )
    commands.append(
        "minimap2 -cx asm5 --secondary=yes consensus_assembly.fasta consensus_assembly.fasta"
    )
    if not wbam.exists():
        p1 = subprocess.Popen(
            [
                tool("minimap2"),
                "-t",
                "12",
                "-ax",
                cfg["reads"]["alignment_preset"],
                mmi,
                paths["reads"],
            ],
            stdout=subprocess.PIPE,
        )
        subprocess.run(
            [tool("samtools"), "sort", "-@", "8", "-m", "1G", "-o", wbam, "-"],
            stdin=p1.stdout,
            check=True,
        )
        p1.stdout.close()
        assert p1.wait() == 0
    if not Path(str(wbam) + ".bai").exists():
        subprocess.run([tool("samtools"), "index", "-@", "8", wbam], check=True)
    commands += [
        "minimap2 -t 12 -ax map-ont source.mmi BC04.trimmed.fastq.gz | "
        "samtools sort -@ 8 -m 1G -o source_reads.bam -",
        "samtools index -@ 8 source_reads.bam",
    ]
    contigs, meta, rels, merges, muts, intervals, cases, junctions, segments = build_sequences(
        seq, cfg
    )
    # deterministic FASTAs
    for sample in ("S01", "S02", "S03"):
        rec = sorted([(k, v) for k, (s, v) in contigs.items() if s == sample])
        write_gzip(out / "contigs" / f"{sample}.contigs.fasta.gz", format_fasta(rec))
    # Select source reads at true junction coordinates and representative loci via indexed BAM.
    import pysam

    chosen = {}
    locus_reads = {}
    bam = pysam.AlignmentFile(wbam, "rb")
    ref = bam.references[0]
    n = len(seq)
    native = [j for j in junctions if j["junction_is_true"] == "true"]
    for j in native:
        c = int(j["source_coordinate"])
        regs = [(max(0, c - 5000), min(n, c + 5000))] if c else [(0, 5000), (n - 5000, n)]
        candidates = {}
        for s, e in regs:
            for r in bam.fetch(ref, s, e):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                name = r.query_name
                span = (
                    (r.reference_start <= max(0, c - 500) and r.reference_end >= min(n, c + 500))
                    if c
                    else (r.reference_start < 5000 or r.reference_end > n - 5000)
                )
                candidates[name] = (span, r)
        ordered = sorted(
            candidates, key=lambda x: (not candidates[x][0], stable_hash(x, cfg["seed"]))
        )[: cfg["reads"]["maximum_reads_per_locus"]]
        sp = sum(candidates[x][0] for x in ordered)
        j["selected_spanning_reads"] = sp
        j["selected_nonspanning_reads"] = len(ordered) - sp
        locus_reads[j["case_id"]] = ordered
        for name in ordered:
            r = candidates[name][1]
            if r.query_sequence and r.qual:
                chosen.setdefault(name, (r.query_sequence, r.qual))
    bam.close()
    samples_reads = {s: [] for s in ("S01", "S02", "S03")}
    for name, (s, q) in chosen.items():
        samples_reads[f"S0{partition(name, 3, cfg['seed']) + 1}"].append((name, s, q))
    for sample, records in samples_reads.items():
        records.sort(key=lambda x: stable_hash(x[0], cfg["seed"]))
        text = "".join(f"@{n}\n{s}\n+\n{q}\n" for n, s, q in records)
        write_gzip(out / "reads" / f"{sample}.targeted.fastq.gz", text)
    # sample BAMs, references exactly sample FASTA
    for sample in ("S01", "S02", "S03"):
        fq = out / "reads" / f"{sample}.targeted.fastq.gz"
        fa = out / "contigs" / f"{sample}.contigs.fasta.gz"
        ob = out / "alignments" / f"{sample}.contigs.bam"
        p = subprocess.Popen(
            [tool("minimap2"), "-t", "4", "-ax", "map-ont", fa, fq], stdout=subprocess.PIPE
        )
        subprocess.run(
            [tool("samtools"), "sort", "-@", "2", "-o", ob, "-"], stdin=p.stdout, check=True
        )
        p.stdout.close()
        assert p.wait() == 0
        subprocess.run([tool("samtools"), "index", ob], check=True)
        commands.append(
            f"minimap2 -ax map-ont contigs/{sample}.contigs.fasta.gz "
            f"reads/{sample}.targeted.fastq.gz | samtools sort -o "
            f"alignments/{sample}.contigs.bam -"
        )
    # all-vs-all
    combined = work / "all_contigs.fasta"
    combined.write_text(format_fasta(sorted((k, v) for k, (s, v) in contigs.items())))
    for preset in ("asm5", "asm20"):
        raw = work / f"all_vs_all.{preset}.paf"
        with raw.open("wb") as f:
            subprocess.run(
                [
                    tool("minimap2"),
                    "-cx",
                    preset,
                    "--secondary=yes",
                    "-p",
                    "0",
                    "-N",
                    "100",
                    combined,
                    combined,
                ],
                stdout=f,
                check=True,
            )
        write_gzip(out / "alignments" / f"all_vs_all.{preset}.paf.gz", raw.read_text())
        commands.append(
            f"minimap2 -cx {preset} --secondary=yes -p 0 -N 100 all_contigs.fasta all_contigs.fasta"
        )
    # Truth and metadata
    rf = [
        "query_id",
        "target_id",
        "expected_relationship",
        "expected_orientation",
        "expected_status",
        "merge_allowed",
        "expected_overlap_length",
        "expected_identity",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "ambiguity_group",
        "case_id",
        "reason",
    ]
    cf = [
        "contig_id",
        "sample",
        "case_id",
        "source_start",
        "source_end",
        "strand",
        "circular_wrap",
        "length",
        "mutation_set",
        "construction_type",
        "description",
    ]
    write_tsv(out / "expected" / "expected_relationships.tsv", rf, rels)
    write_tsv(
        out / "expected" / "expected_contigs.tsv", cf, sorted(meta, key=lambda x: x["contig_id"])
    )
    write_tsv(
        out / "expected" / "expected_contig_segments.tsv",
        ["contig_id", "segment_index", "source_start", "source_end", "strand"],
        sorted(segments, key=lambda x: (x["contig_id"], x["segment_index"])),
    )
    write_tsv(
        out / "expected" / "mutations.tsv",
        [
            "mutation_set",
            "contig_id",
            "position",
            "mutation_type",
            "reference",
            "alternate",
            "case_id",
        ],
        sorted(muts, key=lambda x: (x["mutation_set"], x["position"])),
    )
    write_tsv(
        out / "expected" / "selected_reference_intervals.tsv",
        [
            "case_id",
            "contig_id",
            "source_start",
            "source_end",
            "strand",
            "circular_wrap",
            "selection_rationale",
        ],
        sorted(intervals, key=lambda x: (x["case_id"], x["contig_id"])),
    )
    write_tsv(
        out / "expected" / "expected_junctions.tsv",
        [
            "case_id",
            "left_contig",
            "right_contig",
            "junction_is_true",
            "source_coordinate",
            "circular_wrap",
            "selected_spanning_reads",
            "selected_nonspanning_reads",
            "expected_read_evidence",
            "reason",
        ],
        sorted(junctions, key=lambda x: x["case_id"]),
    )
    comps = [
        {
            "component_id": f"C{i + 1:03d}",
            "case_id": c["case_id"],
            "expected_action": "MERGE"
            if c["category"] == "positive"
            else "KEEP_SEPARATE_OR_AMBIGUOUS",
            "reason": c["description"],
        }
        for i, c in enumerate(sorted(cases, key=lambda x: x["case_id"]))
    ]
    write_tsv(
        out / "expected" / "expected_components.tsv",
        ["component_id", "case_id", "expected_action", "reason"],
        comps,
    )
    write_tsv(
        out / "cases" / "cases.tsv",
        ["case_id", "category", "description", "governing_decision"],
        sorted(cases, key=lambda x: x["case_id"]),
    )
    write_gzip(
        out / "expected" / "expected_merged_sequences.fasta.gz",
        format_fasta(sorted((k, v) for k, v in merges.items())),
    )
    (out / "VERSION").write_text("1.0.0\n")
    (out / "manifest.tsv").write_text(
        "sample\tcontigs\tbam\ttechnology\n"
        + "".join(
            f"{s}\tcontigs/{s}.contigs.fasta.gz\talignments/{s}.contigs.bam\tont\n"
            for s in ("S01", "S02", "S03")
        )
    )
    (out / "manifest_no_reads.tsv").write_text(
        "sample\tcontigs\ttechnology\n"
        + "".join(f"{s}\tcontigs/{s}.contigs.fasta.gz\tont\n" for s in ("S01", "S02", "S03"))
    )
    source["natural_repeat_assessment"] = {
        "self_alignment_file": "private source_self.asm5.paf",
        "decision": (
            "synthetic fallback used for controlled repeat ambiguity; "
            "natural repeat not asserted as truth"
        ),
    }
    write_json(out / "metadata" / "source_metadata.json", source)
    (out / "metadata" / "generation_config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=True))
    versions = []
    for name, argv in [
        ("Python", [sys.executable, "--version"]),
        ("minimap2", [tool("minimap2"), "--version"]),
        ("samtools", [tool("samtools"), "--version"]),
        ("seqkit", [tool("seqkit"), "version"]),
        ("pigz", [tool("pigz"), "--version"]),
    ]:
        z = subprocess.run(argv, text=True, capture_output=True)
        versions.append({"tool": name, "version": ((z.stdout or z.stderr).splitlines() or [""])[0]})
    write_tsv(out / "metadata" / "tool_versions.tsv", ["tool", "version"], versions)
    write_tsv(
        out / "metadata" / "commands.tsv",
        ["index", "command"],
        [{"index": i + 1, "command": x} for i, x in enumerate(commands)],
    )
    (out / "README.md").write_text(
        f"# Contigger compact validation dataset\n\n"
        f"Version 1.0.0; seed {cfg['seed']}. This benchmark tests conservative exact, "
        "containment, overlap, threshold, circular, mutation, negative and "
        "repeat-ambiguity decisions under the rule **a missed merge is preferable to a "
        "false merge**.\n\nIt contains only benchmark-derived reference slices and targeted "
        "ONT reads. It excludes the complete reads, assembly, GFA, source BAM and indexes. "
        "Coordinates are zero-based, half-open; strand and circular wrapping are explicit in "
        "`expected/`. Use `manifest.tsv` with reads/BAMs or `manifest_no_reads.tsv` for "
        "contig-only runs. PAFs are realistic unfiltered asm5/asm20 all-vs-all classifier "
        "inputs. Source-contig BAMs support pileup and extraction; a newly merged sequence "
        "must be remapped before direct junction validation.\n\nAfter unpacking beside "
        "`benchmarking/`, run `python benchmarking/validate_benchmark.py test_data`.\n"
    )
    # Checksums exclude checksum file itself.
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "checksums.sha256")
    (out / "checksums.sha256").write_text(
        "".join(f"{sha256(p)}  {p.relative_to(out).as_posix()}\n" for p in files)
    )
    retained = sum(len(s) for rs in samples_reads.values() for _, s, _ in rs)
    counts = collections.Counter(r["expected_relationship"] for r in rels)
    bench = {
        "benchmark_version": "1.0.0",
        "samples": 3,
        "contigs": len(contigs),
        "cases": len({c["case_id"] for c in cases}),
        "reads_retained": len(chosen),
        "retained_bases": retained,
        "expected_relationship_types": len(counts),
        "relationship_counts": dict(sorted(counts.items())),
        "configured_thresholds": cfg["thresholds"],
        "source_checksums": source["source_sha256"],
        "construction_seed": cfg["seed"],
        "directory_size_bytes": sum(p.stat().st_size for p in out.rglob("*") if p.is_file()),
        "runtime_metadata_excluded_from_determinism": [],
    }
    write_json(out / "metadata" / "benchmark.json", bench)
    # refresh checksums after benchmark.json
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "checksums.sha256")
    (out / "checksums.sha256").write_text(
        "".join(f"{sha256(p)}  {p.relative_to(out).as_posix()}\n" for p in files)
    )
    if not a.keep_work:
        shutil.rmtree(work)
    print(json.dumps(bench, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
