#!/usr/bin/env python3
"""
map_repetitive_elements_pe.py

Run bowtie2 on PE FASTQ pair, parse output in streaming mode, assign each
read pair to its best-scoring repeat-element family, and write a SAM-like
output.

Usage:
    python map_repetitive_elements_pe.py <r1_fastq> <r2_fastq> \
        <bowtie2_db_prefix> <output_rep_sam> <file_list>

Bowtie2 flags (identical to Perl):
    stdbuf -oL bowtie2 -q --sensitive -a -p 3 --no-mixed --no-discordant
        --reorder -x <db> -1 <r1> -2 <r2> 2> <output>.bowtieout

Note: The Perl script does NOT include --no-discordant in the actual command
it runs (the comment-only line has it; the live line does not). We match the
live command exactly.

rRNA special handling and multi-family logic are identical to the SE version.
The key difference: PE uses R1+R2 pairs; alignment score = sum of R1+R2 AS
scores; UMI = first colon-delimited token of the read name.
"""

import sys
import os
import subprocess
import re

RRNA_EXTRA = {
    "RNA28S": "RNA45S",
    "RNA18S": "RNA45S",
    "RNA5-8S": "RNA45S",
    "antisense_RNA28S": "antisense_RNA45S",
    "antisense_RNA18S": "antisense_RNA45S",
    "antisense_RNA5-8S": "antisense_RNA45S",
}
RRNA_EXTRA_REV = {
    "RNA45S": {"RNA28S", "RNA18S", "RNA5-8S"},
    "antisense_RNA45S": {"antisense_RNA28S", "antisense_RNA18S", "antisense_RNA5-8S"},
}

PRINT_BATCH = 10000


def read_filelists(filepath):
    enst2gene = {}
    convert_enst2type = {}
    priority_n = 0
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            parts = line.split("\t")
            allenst = parts[0] if len(parts) > 0 else ""
            type_label = parts[3] if len(parts) > 3 else ""
            gid = parts[2] if len(parts) > 2 else ""
            if not allenst:
                print(f"error missing enst {line} {filepath}", file=sys.stderr)
                continue
            type_label = type_label.rstrip("_")
            for enst in allenst.split("|"):
                enst2gene[enst] = gid
                convert_enst2type[enst] = f"{type_label}:{priority_n}"
                priority_n += 1
    return enst2gene, convert_enst2type


def get_as_score(flags_str):
    m = re.search(r"AS:i:(\S+?)(\s|$)", flags_str)
    if m:
        return int(m.group(1))
    return None


def print_output(read_hash, samout, multimapping_hash):
    for read_name in read_hash:
        rd = read_hash[read_name]
        ensttype_array = sorted(rd["flags"].keys())
        ensttype = ensttype_array[0]

        masterenst_array = [rd["master_enst"][t] for t in ensttype_array]
        ensttype_join = "|".join(ensttype_array)
        masterenst_join = "|".join(masterenst_array)

        if len(rd["flags"]) == 1:
            r1_fields = rd["R1"][ensttype].split("\t")
            r2_fields = rd["R2"][ensttype].split("\t")
            r1_fields[2] = ensttype_join + "||" + masterenst_join
            r2_fields[2] = ensttype_join + "||" + masterenst_join
            r1_line = "\t".join(r1_fields)
            r2_line = "\t".join(r2_fields)
            zz_tag = "|".join(rd["mult_ensts"][ensttype])
            samout.write(r1_line + "\tZZ:Z:" + zz_tag + "\n")
            samout.write(r2_line + "\tZZ:Z:" + zz_tag + "\n")
        else:
            all_mult_ensts = []
            for key in ensttype_array:
                all_mult_ensts.append("|".join(rd["mult_ensts"][key]))
            final_mult_ensts = "|".join(all_mult_ensts)

            r1_fields = rd["R1"][ensttype].split("\t")
            r2_fields = rd["R2"][ensttype].split("\t")
            r1_fields[2] = ensttype_join + "||" + masterenst_join
            r2_fields[2] = ensttype_join + "||" + masterenst_join
            r1_line = "\t".join(r1_fields)
            r2_line = "\t".join(r2_fields)
            samout.write(r1_line + "\tZZ:Z:" + final_mult_ensts + "\n")
            samout.write(r2_line + "\tZZ:Z:" + final_mult_ensts + "\n")

            multimapping_type = "|".join(rd["flags"].keys())
            multimapping_hash[multimapping_type] = (
                multimapping_hash.get(multimapping_type, 0) + 1
            )


def process_alignment(
    r1, r2, r1name, frag_strand,
    paired_mismatch_score, mapped_enst, mapped_enst_full,
    ensttype, enstpriority,
    read_hash, prev_r1name_ref, read_counter_ref,
    samout, multimapping_hash,
):
    if r1name == prev_r1name_ref[0]:
        pass
    else:
        read_counter_ref[0] += 1
        if read_counter_ref[0] > PRINT_BATCH:
            print_output(read_hash, samout, multimapping_hash)
            read_hash.clear()
            read_counter_ref[0] = 0
        prev_r1name_ref[0] = r1name

    if r1name not in read_hash:
        read_hash[r1name] = {
            "R1": {ensttype: r1},
            "R2": {ensttype: r2},
            "flags": {ensttype: enstpriority},
            "quality": paired_mismatch_score,
            "mult_ensts": {ensttype: [mapped_enst_full]},
            "enst": {mapped_enst: mapped_enst_full},
            "master_enst": {ensttype: mapped_enst_full},
        }
    else:
        rd = read_hash[r1name]
        if paired_mismatch_score < rd["quality"]:
            pass  # worse — skip
        elif paired_mismatch_score > rd["quality"]:
            # better — replace everything
            read_hash[r1name] = {
                "R1": {ensttype: r1},
                "R2": {ensttype: r2},
                "flags": {ensttype: enstpriority},
                "quality": paired_mismatch_score,
                "mult_ensts": {ensttype: [mapped_enst_full]},
                "enst": {mapped_enst: mapped_enst_full},
                "master_enst": {ensttype: mapped_enst_full},
            }
        else:
            # equal quality
            if ensttype in rd["flags"]:
                if mapped_enst in rd["enst"]:
                    prev_full = rd["enst"][mapped_enst]
                    if prev_full + "_withgenomeflank" == mapped_enst_full or \
                       prev_full + "_spliced" == mapped_enst_full:
                        pass  # keep old
                    elif rd["enst"][mapped_enst] == mapped_enst_full + "_withgenomeflank" or \
                         rd["enst"][mapped_enst] == mapped_enst_full + "_spliced":
                        rd["R1"][ensttype] = r1
                        rd["R2"][ensttype] = r2
                        rd["flags"][ensttype] = enstpriority
                        for i, v in enumerate(rd["mult_ensts"][ensttype]):
                            if v == prev_full:
                                rd["mult_ensts"][ensttype][i] = mapped_enst_full
                        rd["enst"][mapped_enst] = mapped_enst_full
                        rd["master_enst"][ensttype] = mapped_enst_full
                    else:
                        for i, v in enumerate(rd["mult_ensts"][ensttype]):
                            if v == rd["enst"][mapped_enst]:
                                rd["mult_ensts"][ensttype][i] = mapped_enst_full + "_DOUBLEMAP"
                        rd["master_enst"][ensttype] = mapped_enst_full + "_DOUBLEMAP"
                elif enstpriority < rd["flags"][ensttype]:
                    rd["R1"][ensttype] = r1
                    rd["R2"][ensttype] = r2
                    rd["flags"][ensttype] = enstpriority
                    rd["mult_ensts"][ensttype].insert(0, mapped_enst_full)
                    rd["enst"][mapped_enst] = mapped_enst_full
                    rd["master_enst"][ensttype] = mapped_enst_full
                else:
                    rd["mult_ensts"][ensttype].append(mapped_enst_full)

            elif ensttype in RRNA_EXTRA and RRNA_EXTRA[ensttype] in rd["R1"]:
                old_label = RRNA_EXTRA[ensttype]
                del rd["R1"][old_label]
                del rd["R2"][old_label]
                del rd["flags"][old_label]
                del rd["master_enst"][old_label]
                del rd["mult_ensts"][old_label]
                rd["enst"].pop("NR_046235.1", None)

                rd["R1"][ensttype] = r1
                rd["R2"][ensttype] = r2
                rd["flags"][ensttype] = enstpriority
                rd["quality"] = paired_mismatch_score
                rd["enst"][mapped_enst] = mapped_enst_full
                rd["master_enst"][ensttype] = mapped_enst_full
                rd["mult_ensts"][ensttype] = [mapped_enst_full]

            elif ensttype in RRNA_EXTRA_REV:
                rrrna_flag = any(
                    elem in rd["R1"] for elem in RRNA_EXTRA_REV[ensttype]
                )
                if not rrrna_flag:
                    rd["R1"][ensttype] = r1
                    rd["R2"][ensttype] = r2
                    rd["flags"][ensttype] = enstpriority
                    rd["quality"] = paired_mismatch_score
                    rd["enst"][mapped_enst] = mapped_enst_full
                    rd["mult_ensts"].setdefault(ensttype, []).append(mapped_enst_full)
                    rd["master_enst"][ensttype] = mapped_enst_full
            else:
                rd["R1"][ensttype] = r1
                rd["R2"][ensttype] = r2
                rd["flags"][ensttype] = enstpriority
                rd["quality"] = paired_mismatch_score
                rd["enst"][mapped_enst] = mapped_enst_full
                rd["mult_ensts"].setdefault(ensttype, []).append(mapped_enst_full)
                rd["master_enst"][ensttype] = mapped_enst_full


def check_bowtie_exit(returncode, bowtie_out):
    """Abort if bowtie2 failed.

    Without this the pipe reader treats an immediate EOF as a successful run
    with zero alignments, so a missing bowtie2 or an unreadable index yields an
    empty rep.sam, a '.done' file and exit 0 -- surfacing much later as
    '#READINFO RepFamilyReads 0 0' instead of an error.
    """
    if returncode == 0:
        return
    try:
        with open(bowtie_out) as fh:
            detail = fh.read().strip()
    except OSError:
        detail = "(no stderr captured)"
    sys.exit(
        f"bowtie2 failed with exit code {returncode}; see {bowtie_out}\n{detail}"
    )


def main():
    if len(sys.argv) != 6:
        print(
            "Usage: map_repetitive_elements_pe.py <r1_fastq> <r2_fastq> "
            "<bowtie2_db_prefix> <output_rep_sam> <file_list>",
            file=sys.stderr,
        )
        sys.exit(1)

    fastq_file1 = sys.argv[1]
    fastq_file2 = sys.argv[2]
    bowtie_db = sys.argv[3]
    output = sys.argv[4]
    filelist_file = sys.argv[5]

    enst2gene, convert_enst2type = read_filelists(filelist_file)

    bowtie_out = output + ".bowtieout"
    multimapping_out = output + ".multimapping_deleted"
    done_file = output + ".done"

    bowtie_cmd = (
        f"stdbuf -oL bowtie2 -q --sensitive -a -p 3 --no-mixed --reorder "
        f"-x {bowtie_db} -1 {fastq_file1} -2 {fastq_file2} 2> {bowtie_out}"
    )
    print(f"command {bowtie_cmd}", file=sys.stderr)

    bowtie_err = open(bowtie_out, "w")
    proc = subprocess.Popen(
        ["stdbuf", "-oL", "bowtie2", "-q", "--sensitive", "-a",
         "-p", "3", "--no-mixed", "--reorder",
         "-x", bowtie_db, "-1", fastq_file1, "-2", fastq_file2],
        stdout=subprocess.PIPE,
        stderr=bowtie_err,
        text=True,
    )

    read_hash = {}
    multimapping_hash = {}
    prev_r1name_ref = [""]
    read_counter_ref = [0]

    with open(output, "w") as samout:
        while True:
            r1 = proc.stdout.readline()
            if not r1:
                break
            r1 = r1.rstrip("\n")

            if r1.startswith("@"):
                samout.write(r1 + "\n")
                continue

            r2 = proc.stdout.readline()
            if not r2:
                print(f"PE: missing R2 for R1: {r1}", file=sys.stderr)
                break
            r2 = r2.rstrip("\n")

            tmp_r1 = r1.split("\t")
            tmp_r2 = r2.split("\t")

            r1name = tmp_r1[0].split()[0]
            r2name = tmp_r2[0].split()[0]

            if r1name != r2name:
                print(f"paired end mismatch error: r1 {tmp_r1[0]} r2 {tmp_r2[0]}", file=sys.stderr)

            r1sam_flag = int(tmp_r1[1])

            if r1sam_flag == 77 or r1sam_flag == 141:
                continue

            if r1sam_flag in (99, 355):
                frag_strand = "-"
            elif r1sam_flag in (83, 339):
                frag_strand = "+"
            elif r1sam_flag in (147, 403):
                frag_strand = "-"
            elif r1sam_flag in (163, 419):
                frag_strand = "+"
            else:
                continue

            flags_r1 = "\t".join(tmp_r1[11:]) + " "
            flags_r2 = "\t".join(tmp_r2[11:]) + " "
            score_r1 = get_as_score(flags_r1)
            score_r2 = get_as_score(flags_r2)
            score_r1 = score_r1 if score_r1 is not None else 0
            score_r2 = score_r2 if score_r2 is not None else 0
            paired_mismatch_score = score_r1 + score_r2

            mapped_enst_full = tmp_r1[2]
            mapped_enst = mapped_enst_full
            if "_spliced" in mapped_enst:
                mapped_enst = re.sub(r"_spliced.*$", "", mapped_enst)
            if "_withgenomeflank" in mapped_enst:
                mapped_enst = re.sub(r"_withgenomeflank.*$", "", mapped_enst)

            if mapped_enst not in convert_enst2type:
                print(f"enst2type is missing for {mapped_enst} {r1}", file=sys.stderr)
                continue
            ensttype_raw, enstpriority_s = convert_enst2type[mapped_enst].split(":")
            enstpriority = int(enstpriority_s)
            ensttype = ensttype_raw

            if frag_strand == "-":
                ensttype = "antisense_" + ensttype
                mapped_enst_full = "antisense_" + mapped_enst_full

            process_alignment(
                r1, r2, r1name, frag_strand,
                paired_mismatch_score, mapped_enst, mapped_enst_full,
                ensttype, enstpriority,
                read_hash, prev_r1name_ref, read_counter_ref,
                samout, multimapping_hash,
            )

        print_output(read_hash, samout, multimapping_hash)

    returncode = proc.wait()
    bowtie_err.close()
    check_bowtie_exit(returncode, bowtie_out)

    with open(multimapping_out, "w") as mmout:
        for key in multimapping_hash:
            mmout.write(f"{key}\t{multimapping_hash[key]}\n")

    with open(done_file, "w") as donef:
        donef.write("jobs done\n")


if __name__ == "__main__":
    main()
