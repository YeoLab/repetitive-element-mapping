#!/usr/bin/env python3
"""
downsample_fastq_bam.py

Produce a downsampled FASTQ + BAM dataset where every UMI 2-nt prefix
has at least --min_per_prefix reads. All reads in the output BAM exist
in the output FASTQ (required by the pipeline).

UMI location:
  SE: after the last '_' in the read name   e.g. K00180:..._CGCCTTGCCG -> prefix CG
  PE: before the first ':' in the read name e.g. AGAAA:SN1001:...      -> prefix AG

Usage (SE):
  python downsample_fastq_bam.py --mode SE \\
    --fastq_r1 in.r1.fq.gz --bam in.bam \\
    --out_r1 out.r1.ds.fq.gz --out_bam out.rmRep.ds.bam

Usage (PE):
  python downsample_fastq_bam.py --mode PE \\
    --fastq_r1 in.r1.fq.gz --fastq_r2 in.r2.fq.gz --bam in.bam \\
    --out_r1 out.r1.ds.fq.gz --out_r2 out.r2.ds.fq.gz --out_bam out.rmRep.ds.bam
"""

import argparse
import gzip
import os
import random
import sys
from collections import defaultdict

import pysam

PREFIXES = [
    "AA", "AC", "AG", "AT", "AN",
    "CA", "CC", "CG", "CT", "CN",
    "GA", "GC", "GG", "GT", "GN",
    "TA", "TC", "TG", "TT", "TN",
    "NA", "NC", "NG", "NT", "NN",
]
VALID_BASES = set("ACGTN")


def _norm_prefix(s):
    p0 = s[0] if len(s) > 0 and s[0] in VALID_BASES else "N"
    p1 = s[1] if len(s) > 1 and s[1] in VALID_BASES else "N"
    return p0 + p1


def get_prefix_se(read_name):
    umi = read_name.split("_")[-1]
    return _norm_prefix(umi)


def get_prefix_pe(read_name):
    umi = read_name.split(":")[0]
    return _norm_prefix(umi)


def _bare(read_name):
    return read_name.split("/")[0].split(" ")[0]


def stream_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            name = fh.readline()
            if not name:
                break
            seq  = fh.readline().rstrip()
            plus = fh.readline().rstrip()
            qual = fh.readline().rstrip()
            yield name.rstrip(), seq, plus, qual


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["SE", "PE"], required=True)
    parser.add_argument("--fastq_r1", required=True)
    parser.add_argument("--fastq_r2")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--out_r1", required=True)
    parser.add_argument("--out_r2")
    parser.add_argument("--out_bam", required=True)
    parser.add_argument("--min_per_prefix", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.mode == "PE" and not args.fastq_r2:
        parser.error("--fastq_r2 is required for PE mode")
    if args.mode == "PE" and not args.out_r2:
        parser.error("--out_r2 is required for PE mode")

    random.seed(args.seed)
    get_prefix = get_prefix_se if args.mode == "SE" else get_prefix_pe

    # --- Pass 1: bucket all R1 reads by UMI prefix ---
    buckets = defaultdict(list)   # prefix -> [bare_name, ...]
    r1_recs = {}                   # bare_name -> (name, seq, plus, qual)

    print("Reading R1 FASTQ...", file=sys.stderr)
    for name, seq, plus, qual in stream_fastq(args.fastq_r1):
        bare = _bare(name[1:])    # strip leading @
        prefix = get_prefix(bare)
        buckets[prefix].append(bare)
        r1_recs[bare] = (name, seq, plus, qual)

    # --- Select reads ---
    selected = set()
    for p in PREFIXES:
        pool = buckets.get(p, [])
        if len(pool) == 0:
            print(f"WARNING: prefix {p!r} has 0 reads in R1 FASTQ", file=sys.stderr)
            continue
        n = min(len(pool), args.min_per_prefix)
        selected.update(random.sample(pool, n))

    print(f"Selected {len(selected)} reads total", file=sys.stderr)
    for p in PREFIXES:
        n = sum(1 for nm in buckets.get(p, []) if nm in selected)
        print(f"  {p}: {n}", file=sys.stderr)

    # --- Write R1 ---
    os.makedirs(os.path.dirname(os.path.abspath(args.out_r1)), exist_ok=True)
    print("Writing R1...", file=sys.stderr)
    with gzip.open(args.out_r1, "wt") as fh:
        for bare, (name, seq, plus, qual) in r1_recs.items():
            if bare in selected:
                fh.write(f"{name}\n{seq}\n{plus}\n{qual}\n")

    # --- Write R2 (PE) ---
    if args.mode == "PE":
        print("Writing R2...", file=sys.stderr)
        with gzip.open(args.out_r2, "wt") as fh:
            for name, seq, plus, qual in stream_fastq(args.fastq_r2):
                bare = _bare(name[1:])
                if bare in selected:
                    fh.write(f"{name}\n{seq}\n{plus}\n{qual}\n")

    # --- Filter BAM ---
    print("Filtering BAM...", file=sys.stderr)
    tmp_bam = args.out_bam + ".unsorted.bam"
    with pysam.AlignmentFile(args.bam, "rb") as bam_in:
        with pysam.AlignmentFile(tmp_bam, "wb", header=bam_in.header) as bam_out:
            kept = 0
            for read in bam_in.fetch(until_eof=True):
                if read.query_name in selected:
                    bam_out.write(read)
                    kept += 1
    print(f"  {kept} BAM alignments kept", file=sys.stderr)

    pysam.sort("-o", args.out_bam, tmp_bam)
    os.remove(tmp_bam)
    pysam.index(args.out_bam)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
