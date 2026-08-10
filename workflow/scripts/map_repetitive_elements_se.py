#!/usr/bin/env python3
"""
map_repetitive_elements_se.py

Run bowtie2 on SE FASTQ, parse output in streaming mode, assign each read to
its best-scoring repeat-element family, and write a SAM-like output.

Usage:
    python map_repetitive_elements_se.py <r1_fastq> <bowtie2_db_prefix> \
        <output_rep_sam> <file_list>

Bowtie2 flags (identical to Perl):
    stdbuf -oL bowtie2 -q --sensitive -a -p 3 --no-mixed --reorder
        -x <db> -U <r1> 2> <output>.bowtieout

rRNA special handling:
    RNA28S, RNA18S, RNA5-8S  -> RNA45S  (forward)
    antisense_ variants      -> antisense_RNA45S

Multi-family reads: when a read maps equally well to multiple families, all
family names are joined with "|" in the REFERENCE field.
"""

import sys
import os
import subprocess
import re

# ---------------------------------------------------------------------------
# rRNA remapping tables (identical to Perl)
# ---------------------------------------------------------------------------
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
    """Return (enst2gene, convert_enst2type) dicts from a 5-column TSV."""
    enst2gene = {}
    convert_enst2type = {}
    priority_n = 0
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            parts = line.split("\t")
            allenst = parts[0] if len(parts) > 0 else ""
            allensg = parts[1] if len(parts) > 1 else ""
            gid = parts[2] if len(parts) > 2 else ""
            type_label = parts[3] if len(parts) > 3 else ""
            # typefile = parts[4]  # not used here
            type_label = type_label.rstrip("_")
            if not allenst:
                print(f"error missing enst {line} {filepath}", file=sys.stderr)
                continue
            for enst in allenst.split("|"):
                enst2gene[enst] = gid
                convert_enst2type[enst] = f"{type_label}:{priority_n}"
                priority_n += 1
    return enst2gene, convert_enst2type


def get_as_score(flags_str):
    """Extract AS:i:<N> from the optional-field string. Returns int or None."""
    m = re.search(r"AS:i:(\S+?)(\s|$)", flags_str)
    if m:
        return int(m.group(1))
    return None


def print_output(read_hash, samout, multimapping_hash):
    """Flush current read_hash to samout, update multimapping_hash in place."""
    for read_name in read_hash:
        rd = read_hash[read_name]
        # Sort family types lexicographically — matches Perl `sort {$a cmp $b}`
        ensttype_array = sorted(rd["flags"].keys())
        ensttype = ensttype_array[0]

        masterenst_array = [rd["master_enst"][t] for t in ensttype_array]
        ensttype_join = "|".join(ensttype_array)
        masterenst_join = "|".join(masterenst_array)

        if len(rd["flags"]) == 1:
            r1_fields = rd["R1"][ensttype].split("\t")
            r1_fields[2] = ensttype_join + "||" + masterenst_join
            r1_line = "\t".join(r1_fields)
            zz_tag = "|".join(rd["mult_ensts"][ensttype])
            samout.write(r1_line + "\tZZ:Z:" + zz_tag + "\n")
        else:
            all_mult_ensts = []
            for key in ensttype_array:
                all_mult_ensts.append("|".join(rd["mult_ensts"][key]))
            final_mult_ensts = "|".join(all_mult_ensts)

            r1_fields = rd["R1"][ensttype].split("\t")
            r1_fields[2] = ensttype_join + "||" + masterenst_join
            r1_line = "\t".join(r1_fields)
            samout.write(r1_line + "\tZZ:Z:" + final_mult_ensts + "\n")

            multimapping_type = "|".join(rd["flags"].keys())
            multimapping_hash[multimapping_type] = (
                multimapping_hash.get(multimapping_type, 0) + 1
            )


def process_alignment(
    r1, r1name, r1sam_flag, frag_strand,
    flags_str, mismatch_score, mapped_enst, mapped_enst_full,
    ensttype, enstpriority,
    read_hash, prev_r1name_ref, read_counter_ref,
    samout, multimapping_hash,
):
    """
    Core logic: update read_hash for the current alignment.
    prev_r1name_ref and read_counter_ref are single-element lists (mutable refs).
    """
    r1name_cur = r1name
    paired_mismatch_score = mismatch_score

    # Batch flush
    if r1name_cur == prev_r1name_ref[0]:
        pass
    else:
        read_counter_ref[0] += 1
        if read_counter_ref[0] > PRINT_BATCH:
            print_output(read_hash, samout, multimapping_hash)
            read_hash.clear()
            read_counter_ref[0] = 0
        prev_r1name_ref[0] = r1name_cur

    if r1name_cur not in read_hash:
        read_hash[r1name_cur] = {
            "R1": {ensttype: r1},
            "flags": {ensttype: enstpriority},
            "quality": paired_mismatch_score,
            "mult_ensts": {ensttype: [mapped_enst_full]},
            "enst": {mapped_enst: mapped_enst_full},
            "master_enst": {ensttype: mapped_enst_full},
        }
    else:
        rd = read_hash[r1name_cur]
        if paired_mismatch_score < rd["quality"]:
            # Worse match — skip
            pass
        elif paired_mismatch_score > rd["quality"]:
            # Better match — replace everything
            read_hash[r1name_cur] = {
                "R1": {ensttype: r1},
                "flags": {ensttype: enstpriority},
                "quality": paired_mismatch_score,
                "mult_ensts": {ensttype: [mapped_enst_full]},
                "enst": {mapped_enst: mapped_enst_full},
                "master_enst": {ensttype: mapped_enst_full},
            }
        else:
            # Equal quality
            if ensttype in rd["flags"]:
                # Same family exists already
                if mapped_enst in rd["enst"]:
                    prev_full = rd["enst"][mapped_enst]
                    if prev_full + "_withgenomeflank" == mapped_enst_full or \
                       prev_full + "_spliced" == mapped_enst_full:
                        # Keep old (shorter) — skip new
                        pass
                    elif rd["enst"][mapped_enst] == mapped_enst_full + "_withgenomeflank" or \
                         rd["enst"][mapped_enst] == mapped_enst_full + "_spliced":
                        # Replace old longer with new shorter
                        rd["R1"][ensttype] = r1
                        rd["flags"][ensttype] = enstpriority
                        for i, v in enumerate(rd["mult_ensts"][ensttype]):
                            if v == prev_full:
                                rd["mult_ensts"][ensttype][i] = mapped_enst_full
                        rd["enst"][mapped_enst] = mapped_enst_full
                        rd["master_enst"][ensttype] = mapped_enst_full
                    else:
                        # Double map — flag
                        for i, v in enumerate(rd["mult_ensts"][ensttype]):
                            if v == rd["enst"][mapped_enst]:
                                rd["mult_ensts"][ensttype][i] = mapped_enst_full + "_DOUBLEMAP"
                        rd["master_enst"][ensttype] = mapped_enst_full + "_DOUBLEMAP"
                elif enstpriority < rd["flags"][ensttype]:
                    # Better priority within same family
                    rd["R1"][ensttype] = r1
                    rd["flags"][ensttype] = enstpriority
                    rd["mult_ensts"][ensttype].insert(0, mapped_enst_full)
                    rd["enst"][mapped_enst] = mapped_enst_full
                    rd["master_enst"][ensttype] = mapped_enst_full
                else:
                    # Same or worse priority — keep new enst_full
                    rd["mult_ensts"][ensttype].append(mapped_enst_full)

            elif ensttype in RRNA_EXTRA and RRNA_EXTRA[ensttype] in rd["R1"]:
                # New is specific rRNA (RNA28S etc); old was generic rRNA45S — replace
                old_label = RRNA_EXTRA[ensttype]
                del rd["R1"][old_label]
                del rd["flags"][old_label]
                del rd["master_enst"][old_label]
                del rd["mult_ensts"][old_label]
                rd["enst"].pop("NR_046235.1", None)

                rd["R1"][ensttype] = r1
                rd["flags"][ensttype] = enstpriority
                rd["quality"] = paired_mismatch_score
                rd["enst"][mapped_enst] = mapped_enst_full
                rd["master_enst"][ensttype] = mapped_enst_full
                rd["mult_ensts"][ensttype] = [mapped_enst_full]

            elif ensttype in RRNA_EXTRA_REV:
                # New is generic rRNA45S
                rrrna_flag = any(
                    elem in rd["R1"] for elem in RRNA_EXTRA_REV[ensttype]
                )
                if not rrrna_flag:
                    # Old is something non-rRNA — multi-family case
                    rd["R1"][ensttype] = r1
                    rd["flags"][ensttype] = enstpriority
                    rd["quality"] = paired_mismatch_score
                    rd["enst"][mapped_enst] = mapped_enst_full
                    rd["mult_ensts"].setdefault(ensttype, []).append(mapped_enst_full)
                    rd["master_enst"][ensttype] = mapped_enst_full
            else:
                # Maps to a new family — multi-family
                rd["R1"][ensttype] = r1
                rd["flags"][ensttype] = enstpriority
                rd["quality"] = paired_mismatch_score
                rd["enst"][mapped_enst] = mapped_enst_full
                rd["mult_ensts"].setdefault(ensttype, []).append(mapped_enst_full)
                rd["master_enst"][ensttype] = mapped_enst_full


def main():
    if len(sys.argv) != 5:
        print(
            "Usage: map_repetitive_elements_se.py <r1_fastq> <bowtie2_db_prefix> "
            "<output_rep_sam> <file_list>",
            file=sys.stderr,
        )
        sys.exit(1)

    fastq_file1 = sys.argv[1]
    bowtie_db = sys.argv[2]
    output = sys.argv[3]
    filelist_file = sys.argv[4]

    enst2gene, convert_enst2type = read_filelists(filelist_file)

    bowtie_out = output + ".bowtieout"
    multimapping_out = output + ".multimapping_deleted"
    done_file = output + ".done"

    bowtie_cmd = (
        f"stdbuf -oL bowtie2 -q --sensitive -a -p 3 --no-mixed --reorder "
        f"-x {bowtie_db} -U {fastq_file1} 2> {bowtie_out}"
    )
    print(f"command {bowtie_cmd}", file=sys.stderr)

    proc = subprocess.Popen(
        ["stdbuf", "-oL", "bowtie2", "-q", "--sensitive", "-a",
         "-p", "3", "--no-mixed", "--reorder",
         "-x", bowtie_db, "-U", fastq_file1],
        stdout=subprocess.PIPE,
        stderr=open(bowtie_out, "w"),
        text=True,
    )

    read_hash = {}
    multimapping_hash = {}
    prev_r1name_ref = [""]
    read_counter_ref = [0]

    with open(output, "w") as samout:
        for r1 in proc.stdout:
            r1 = r1.rstrip("\n")

            if r1.startswith("@"):
                samout.write(r1 + "\n")
                continue

            tmp_r1 = r1.split("\t")

            # UMI extraction for SE: last underscore-delimited token
            read_name_parts = tmp_r1[0].split("_")
            r1bc = read_name_parts[-1]
            r1name = "_".join(read_name_parts[:-1])

            r1sam_flag = int(tmp_r1[1])
            if r1sam_flag == 4:
                continue

            if r1sam_flag in (16, 272):
                frag_strand = "-"
            elif r1sam_flag in (0, 256):
                frag_strand = "+"
            else:
                continue

            flags_str = "\t".join(tmp_r1[11:])
            mismatch_score = get_as_score(flags_str + " ")
            if mismatch_score is None:
                mismatch_score = 0

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
                r1, r1name, r1sam_flag, frag_strand,
                flags_str, mismatch_score, mapped_enst, mapped_enst_full,
                ensttype, enstpriority,
                read_hash, prev_r1name_ref, read_counter_ref,
                samout, multimapping_hash,
            )

        # Final flush
        print_output(read_hash, samout, multimapping_hash)

    proc.wait()

    with open(multimapping_out, "w") as mmout:
        for key in multimapping_hash:
            mmout.write(f"{key}\t{multimapping_hash[key]}\n")

    with open(done_file, "w") as donef:
        donef.write("jobs done\n")


if __name__ == "__main__":
    main()
