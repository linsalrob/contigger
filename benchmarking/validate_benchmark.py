#!/usr/bin/env python3
import argparse,csv,gzip,hashlib,json,os,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from benchmarklib.fasta import read_fasta
from benchmarklib.utilities import sha256

REQUIRED=["README.md","VERSION","manifest.tsv","manifest_no_reads.tsv","metadata/benchmark.json","metadata/source_metadata.json","metadata/generation_config.yaml","metadata/tool_versions.tsv","metadata/commands.tsv","expected/expected_relationships.tsv","expected/expected_contigs.tsv","expected/expected_junctions.tsv","expected/expected_components.tsv","expected/expected_merged_sequences.fasta.gz","expected/mutations.tsv","expected/selected_reference_intervals.tsv","cases/cases.tsv","checksums.sha256"]+[f"contigs/S0{i}.contigs.fasta.gz" for i in range(1,4)]+[f"reads/S0{i}.targeted.fastq.gz" for i in range(1,4)]+[f"alignments/S0{i}.contigs.bam" for i in range(1,4)]+[f"alignments/S0{i}.contigs.bam.bai" for i in range(1,4)]+["alignments/all_vs_all.asm5.paf.gz","alignments/all_vs_all.asm20.paf.gz"]
def tsv(p):
    with open(p,newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def fail(errors,msg):errors.append(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument("dataset");ap.add_argument("--archive");a=ap.parse_args();root=Path(a.dataset).resolve();errors=[]
    for x in REQUIRED:
        if not (root/x).is_file():fail(errors,f"missing: {x}")
    if errors: print("VALIDATION FAILED\n"+"\n".join(errors),file=sys.stderr);return 1
    forbidden={"BC04.trimmed.fastq.gz","consensus_assembly.fasta","consensus_assembly.gfa","source_reads.bam","source.mmi"}
    for p in root.rglob("*"):
        if p.name in forbidden:fail(errors,f"forbidden source/intermediate copied: {p.relative_to(root)}")
        if p.is_symlink():fail(errors,f"symlink forbidden: {p.relative_to(root)}")
        if p.is_file() and p.stat().st_size>50*1024*1024:fail(errors,f"file exceeds 50 MiB: {p.relative_to(root)}")
    contigs={};sample_ids={}
    for i in range(1,4):
        sample=f"S0{i}"; rec=read_fasta(root/"contigs"/f"{sample}.contigs.fasta.gz");sample_ids[sample]={n:len(s) for n,s in rec}
        for n,s in rec:
            if n in contigs:fail(errors,f"duplicate contig id: {n}")
            contigs[n]=s
            if len(s)>=6_852_769:fail(errors,f"full-length assembly-like sequence: {n}")
    truth=tsv(root/"expected"/"expected_contigs.tsv");truthids=[r["contig_id"] for r in truth]
    if set(truthids)!=set(contigs) or len(truthids)!=len(set(truthids)):fail(errors,"truth contigs do not match globally unique FASTA contigs")
    for r in truth:
        if r["contig_id"] in contigs and int(r["length"])!=len(contigs[r["contig_id"]]):fail(errors,f"length mismatch: {r['contig_id']}")
        if r["circular_wrap"]=="true":
            if r["source_start"] and int(r["source_end"])-int(r["source_start"])!=int(r["length"]):fail(errors,f"invalid circular interval: {r['contig_id']}")
    samtools=os.environ.get("SAMTOOLS","/home/edwa0468/miniconda3/envs/contigger-benchmark/bin/samtools")
    for sample,ids in sample_ids.items():
        bam=root/"alignments"/f"{sample}.contigs.bam"
        z=subprocess.run([samtools,"quickcheck",str(bam)],capture_output=True)
        if z.returncode:fail(errors,f"BAM quickcheck failed: {sample}")
        hdr=subprocess.check_output([samtools,"view","-H",str(bam)],text=True);sq={}
        for line in hdr.splitlines():
            if line.startswith("@HD") and "SO:coordinate" not in line:fail(errors,f"BAM not coordinate sorted: {sample}")
            if line.startswith("@SQ"):
                d=dict(x.split(":",1) for x in line.split("\t")[1:]);sq[d["SN"]]=int(d["LN"])
        if sq!=ids:fail(errors,f"BAM dictionary differs from sample FASTA: {sample}")
    readnames=set()
    for i in range(1,4):
        with gzip.open(root/"reads"/f"S0{i}.targeted.fastq.gz","rt") as f:
            while True:
                h=f.readline()
                if not h:break
                s=f.readline();p=f.readline();q=f.readline()
                if not(s and p and q and h.startswith("@") and p.startswith("+") and len(s.rstrip())==len(q.rstrip())):fail(errors,f"malformed FASTQ S0{i}");break
                n=h[1:].split()[0]
                if n in readnames:fail(errors,f"read occurs in multiple FASTQs: {n}")
                readnames.add(n)
    rel=tsv(root/"expected"/"expected_relationships.tsv");classes={r["expected_relationship"] for r in rel};required_classes={"EXACT_MATCH","QUERY_CONTAINED_IN_TARGET","TARGET_CONTAINED_IN_QUERY","QUERY_SUFFIX_TO_TARGET_PREFIX","TARGET_SUFFIX_TO_QUERY_PREFIX","AMBIGUOUS_OVERLAP","NO_RELATIONSHIP"}
    if not required_classes<=classes:fail(errors,f"missing relationship classes: {sorted(required_classes-classes)}")
    for r in rel:
        if r["query_id"] not in contigs or r["target_id"] not in contigs:fail(errors,f"relationship references absent contig: {r['case_id']}")
    muts=tsv(root/"expected"/"mutations.tsv")
    expected={"identity_9801":199,"identity_9800":200,"identity_9799":201}
    for case,n in expected.items():
        got=sum(r["case_id"]==case and r["mutation_type"]=="substitution" for r in muts)
        if got!=n:fail(errors,f"mutation count {case}: {got} != {n}")
    merged=dict(read_fasta(root/"expected"/"expected_merged_sequences.fasta.gz"))
    for case,s in merged.items():
        rows=[r for r in rel if r["case_id"]==case and r["merge_allowed"]=="true" and r["expected_relationship"]=="QUERY_SUFFIX_TO_TARGET_PREFIX"]
        if rows:
            r=rows[0];ov=int(r["expected_overlap_length"]);q=contigs[r["query_id"]];t=contigs[r["target_id"]]
            if r["expected_orientation"]=="+" and q+t[ov:]!=s:fail(errors,f"merged sequence mismatch: {case}")
    for preset in ("asm5","asm20"):
        last=None
        with gzip.open(root/"alignments"/f"all_vs_all.{preset}.paf.gz","rt") as f:
            for no,line in enumerate(f,1):
                z=line.rstrip().split("\t")
                if len(z)<12:fail(errors,f"invalid {preset} PAF line {no}");break
                try:[int(z[i]) for i in (1,2,3,6,7,8,9,10,11)]
                except ValueError:fail(errors,f"non-numeric {preset} PAF line {no}")
    listed=[]
    for line in (root/"checksums.sha256").read_text().splitlines():
        h,name=line.split("  ",1);listed.append(name)
        if not (root/name).is_file() or sha256(root/name)!=h:fail(errors,f"checksum mismatch: {name}")
    actual=sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and p.name!="checksums.sha256")
    if listed!=actual:fail(errors,"checksums list is incomplete or not deterministically sorted")
    # deterministic ordering of primary truth tables
    if truthids!=sorted(truthids):fail(errors,"expected_contigs.tsv ordering is not deterministic")
    relkeys=[(r["case_id"],r["query_id"],r["target_id"]) for r in rel]
    if relkeys!=sorted(relkeys):fail(errors,"expected_relationships.tsv ordering is not deterministic")
    archive=Path(a.archive) if a.archive else root.parent/(root.name+".tar.gz")
    if archive.exists() and archive.stat().st_size>150*1024*1024:fail(errors,"archive exceeds 150 MiB")
    if errors:
        print("VALIDATION FAILED",file=sys.stderr);print("\n".join(f"- {e}" for e in errors),file=sys.stderr);return 1
    print(f"VALIDATION PASSED: {len(contigs)} contigs, {len(readnames)} unique targeted reads, {len(rel)} relationships")
    return 0
if __name__=="__main__":raise SystemExit(main())
