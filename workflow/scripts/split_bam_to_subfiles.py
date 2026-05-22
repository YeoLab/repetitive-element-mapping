#!/usr/bin/env python3
"""
split_bam_to_subfiles.py

Split a SAM/BAM file into 25 subfiles by the first 2 nucleotides of the UMI
randomer. Output files are written to CWD:
    {prefix}.{basename}.tmp

Usage:
    python split_bam_to_subfiles.py <sam_or_bam_path> <SE_or_PE>

For SE reads: UMI is the last underscore-delimited token in the read name.
For PE reads: UMI is the first colon-delimited token in the read name.

After the script runs, the Snakemake rule renames:
    {prefix}.rep.sam.tmp  -> {prefix}.rep.tmp
    {prefix}.rmrep.sam.tmp -> {prefix}.rmrep.tmp
"""

import sys
import os
import subprocess

def main():
    if len(sys.argv) != 3:
        print("Usage: split_bam_to_subfiles.py <sam_file> <SE_or_PE>", file=sys.stderr)
        sys.exit(1)

    sam_fi = sys.argv[1]
    se_or_pe = sys.argv[2]

    if se_or_pe == "PE":
        print("splitting in paired-end mode", file=sys.stderr)
    elif se_or_pe == "SE":
        print("splitting in single-end mode", file=sys.stderr)
    else:
        print("fatal error - SE or PE not defined", file=sys.stderr)
        sys.exit(1)

    # basename of input file (no directory component)
    sam_fi_short = os.path.basename(sam_fi)

    # Open all 25 output file handles in CWD
    bases = ["A", "C", "G", "T", "N"]
    filehandles = {}
    for b1 in bases:
        for b2 in bases:
            prefix = b1 + b2
            outfi = prefix + "." + sam_fi_short + ".tmp"
            filehandles[prefix] = open(outfi, "w")

    # Open the input — SAM directly, BAM via samtools
    if sam_fi.endswith(".sam"):
        infile = open(sam_fi, "r")
    elif sam_fi.endswith(".bam"):
        proc = subprocess.Popen(
            ["samtools", "view", "-h", sam_fi],
            stdout=subprocess.PIPE,
            text=True,
        )
        infile = proc.stdout
    else:
        print(f"weird - {sam_fi} not either sam or bam file format - exit", file=sys.stderr)
        sys.exit(1)

    if se_or_pe == "PE":
        # Read pairs in lockstep
        while True:
            r1 = infile.readline()
            if not r1:
                break
            r1 = r1.rstrip("\n")

            # Skip header lines
            if r1.startswith("@"):
                continue

            r2 = infile.readline()
            if not r2:
                print(f"PE mode: missing R2 for R1: {r1}", file=sys.stderr)
                break
            r2 = r2.rstrip("\n")

            tmp_r1 = r1.split("\t")
            tmp_r2 = r2.split("\t")

            r1name = tmp_r1[0].split()[0]
            r2name = tmp_r2[0].split()[0]
            r1sam_flag = int(tmp_r1[1])

            if r1name != r2name:
                print(f"paired end mismatch error: {sam_fi} r1 {tmp_r1[0]} r2 {tmp_r2[0]}", file=sys.stderr)

            # Skip unmapped pairs
            if r1sam_flag == 77 or r1sam_flag == 141:
                continue

            # Determine fragment strand; also swap R1/R2 for certain flag values
            if r1sam_flag in (99, 355):
                frag_strand = "-"
            elif r1sam_flag in (83, 339):
                frag_strand = "+"
            elif r1sam_flag in (147, 403):
                frag_strand = "-"
                tmp_r1, tmp_r2 = tmp_r2, tmp_r1
                r1, r2 = r2, r1
            elif r1sam_flag in (163, 419):
                frag_strand = "+"
                tmp_r1, tmp_r2 = tmp_r2, tmp_r1
                r1, r2 = r2, r1
            else:
                # Unknown flag — skip
                continue

            # UMI: first colon-delimited token of read name
            randommer = tmp_r1[0].split(":")[0]
            first2rand = randommer[:2]

            if first2rand in filehandles:
                filehandles[first2rand].write(r1 + "\n" + r2 + "\n")
            else:
                print(f"unexpected UMI prefix {first2rand!r} in read {tmp_r1[0]}", file=sys.stderr)

    else:  # SE
        for r1 in infile:
            r1 = r1.rstrip("\n")

            if r1.startswith("@"):
                continue

            tmp_r1 = r1.split("\t")
            r1sam_flag = int(tmp_r1[1])

            # Skip unmapped
            if r1sam_flag == 4:
                continue

            # Determine strand
            if r1sam_flag in (16, 272):
                frag_strand = "-"
            elif r1sam_flag in (0, 256):
                frag_strand = "+"
            else:
                continue

            # UMI: last underscore-delimited token of read name
            read_name_parts = tmp_r1[0].split("_")
            randommer = read_name_parts[-1]
            first2rand = randommer[:2]

            if first2rand in filehandles:
                filehandles[first2rand].write(r1 + "\n")
            else:
                print(f"unexpected UMI prefix {first2rand!r} in read {tmp_r1[0]}", file=sys.stderr)

    infile.close()
    for fh in filehandles.values():
        fh.close()


if __name__ == "__main__":
    main()
