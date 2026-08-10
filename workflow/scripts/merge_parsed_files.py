#!/usr/bin/env python3
"""
merge_parsed_files.py

Merge multiple .parsed_v2.20201210.txt files into one, summing read counts
per element across all input files.

Usage:
    python merge_parsed_files.py <output_path> <input1> <input2> ...

Input format (per file):
    #READINFO lines (two kinds):
        #READINFO\tAll reads:\t<N>\tPCR duplicates removed:\t...  (legacy format)
        #READINFO\tUsableReads\t<N>
        #READINFO\tGenomicReads\t<N>\t<fraction>
        #READINFO\tRepFamilyReads\t<N>\t<fraction>
    TOTAL lines:
        TOTAL\t<element>\t<readnum>\t<rpm>
    ELEMENT lines:
        ELEMENT\t<ensg_primary>\t<readnum>\t<rpm>\t<enst_all>\t<ensg_all>

Output format:
    #READINFO\tAllReads\t<N>
    #READINFO\tUsableReads\t<N>\t<fraction>
    #READINFO\tGenomicReads\t<N>\t<fraction>
    #READINFO\tRepFamilyReads\t<N>\t<fraction>
    TOTAL lines sorted by readnum descending
    ELEMENT lines sorted by readnum descending
"""

import sys
import re


def parse_file(path, read_sums):
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            if not tmp:
                continue

            if tmp[0] == "#READINFO":
                # Legacy format: matches the long "All reads:" pattern
                m = re.match(
                    r"#READINFO\tAll reads:\t(\d+)\tPCR duplicates removed:\t(\d+)"
                    r"\tUsable Remaining:\t(\d+)\tUsable from genomic mapping:\t(\d+)"
                    r"\tUsable from family mapping:\t(\d+)$",
                    line,
                )
                if m:
                    read_sums["all"] = read_sums.get("all", 0) + int(m.group(1))
                elif len(tmp) >= 3 and tmp[1] == "UsableReads":
                    read_sums["usable"] = read_sums.get("usable", 0) + int(tmp[2])
                elif len(tmp) >= 3 and tmp[1] == "GenomicReads":
                    read_sums["genomic"] = read_sums.get("genomic", 0) + int(tmp[2])
                elif len(tmp) >= 3 and tmp[1] == "RepFamilyReads":
                    read_sums["repfamily"] = read_sums.get("repfamily", 0) + int(tmp[2])
                elif len(tmp) >= 3 and tmp[1] == "AllReads":
                    # New format #READINFO\tAllReads\t<N>
                    read_sums["all"] = read_sums.get("all", 0) + int(tmp[2])
                else:
                    print(f"couldn't parse readinfo line {line}", file=sys.stderr)

            elif tmp[0] == "#READINFO2":
                pass  # ignored

            elif tmp[0] == "TOTAL":
                element = tmp[1]
                readnum = int(tmp[2])
                read_sums.setdefault("total", {})
                read_sums["total"][element] = read_sums["total"].get(element, 0) + readnum

            else:
                # ELEMENT line
                # format: ELEMENT\tensg_primary\treadnum\trpm\tenst_all\tensg_all
                if len(tmp) < 6:
                    print(f"short ELEMENT line: {line}", file=sys.stderr)
                    continue
                _ele_delete, ensg_primary, readnum_s, rpm, enst_all, ensg_all = (
                    tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5]
                )
                readnum = int(readnum_s)
                read_sums.setdefault("element", {})
                el = read_sums["element"]
                if enst_all in el:
                    if el[enst_all]["ensg_all"] != ensg_all:
                        print(
                            f"error - ensg_all mismatch {ensg_all} "
                            f"{el[enst_all]['ensg_all']}",
                            file=sys.stderr,
                        )
                    if el[enst_all]["ensg_primary"] != ensg_primary:
                        print(
                            f"error - ensg_primary mismatch {ensg_primary} "
                            f"{el[enst_all]['ensg_primary']}",
                            file=sys.stderr,
                        )
                else:
                    el[enst_all] = {"ensg_all": ensg_all, "ensg_primary": ensg_primary, "readnum": 0}
                el[enst_all]["readnum"] += readnum


def main():
    if len(sys.argv) < 3:
        print("Usage: merge_parsed_files.py <output> <input1> [input2 ...]", file=sys.stderr)
        sys.exit(1)

    output_fi = sys.argv[1]
    files = sys.argv[2:]

    # For legacy compatibility: print short filename and working dir of first file
    import os
    split_fi1 = files[0].split("/")
    short_fi1 = split_fi1[-1]
    working_dir = "/".join(split_fi1[:-1]) if len(split_fi1) > 1 else "."
    print(f"SHORTFI:{short_fi1}", end="")
    print(f"WORKDIR:{working_dir}", end="")

    read_sums = {}
    for f in files:
        parse_file(f, read_sums)

    usable = read_sums.get("usable", 0)

    with open(output_fi, "w") as out:
        out.write(f"#READINFO\tAllReads\t{read_sums.get('all', 0)}\n")
        if usable > 0:
            out.write(
                f"#READINFO\tUsableReads\t{usable}\t"
                f"{usable / usable:.5f}\n"
            )
            genomic = read_sums.get("genomic", 0)
            repfamily = read_sums.get("repfamily", 0)
            out.write(
                f"#READINFO\tGenomicReads\t{genomic}\t"
                f"{genomic / usable:.5f}\n"
            )
            out.write(
                f"#READINFO\tRepFamilyReads\t{repfamily}\t"
                f"{repfamily / usable:.5f}\n"
            )
        else:
            out.write(f"#READINFO\tUsableReads\t0\t0.0\n")
            out.write(f"#READINFO\tGenomicReads\t0\t0.0\n")
            out.write(f"#READINFO\tRepFamilyReads\t0\t0.0\n")

        # TOTAL lines sorted by readnum descending
        total = read_sums.get("total", {})
        for element in sorted(total, key=lambda k: -total[k]):
            rpm = total[element] / usable if usable > 0 else 0
            out.write(f"TOTAL\t{element}\t{total[element]}\t{rpm:.5f}\n")

        # ELEMENT lines sorted by readnum descending
        element_dict = read_sums.get("element", {})
        for enst_all in sorted(element_dict, key=lambda k: -element_dict[k]["readnum"]):
            el = element_dict[enst_all]
            rpm = el["readnum"] / usable if usable > 0 else 0
            out.write(
                f"ELEMENT\t{el['ensg_primary']}\t{el['readnum']}\t"
                f"{rpm:.5f}\t{enst_all}\t{el['ensg_all']}\n"
            )


if __name__ == "__main__":
    main()
