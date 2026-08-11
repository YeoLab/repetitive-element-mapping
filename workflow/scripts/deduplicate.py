#!/usr/bin/env python3
"""
deduplicate.py

UMI-based deduplication with conflict resolution between repeat-family and
unique-genomic mappings. Translates the Perl script:
  duplicate_removal_inline_paired...20201210_simple.pl

Usage:
    python deduplicate.py <rep_tmp> <rmrep_tmp> <SE_or_PE> \
        <gencode_gtf> <gencode_table_browser> <rep_mask_bed> <file_list>

Outputs (written to CWD, named from basename of <rep_tmp>):
    {base}.combined_w_uniquemap.rmDup.sam
    {base}.combined_w_uniquemap.prermDup.sam
    {base}.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt
    {base}.combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt.done

Key decisions:
- Hash iteration order: all dict/set iterations that affect output use
  sorted() to match the Perl 5.18+ mitigation already present in the Perl.
- Conflict threshold: unique genome wins if genome score > rep score + 24.
- CIGAR parsing: 1-based closed-ended SAM positions converted to 0-based
  half-open intervals matching the Perl logic exactly.
"""

import sys
import os
import re
import math
import subprocess
from collections import defaultdict

GENOME_HASHING_VALUE = 1000
REVSTRAND = {"+": "-", "-": "+"}


# ---------------------------------------------------------------------------
# Reference data loading
# ---------------------------------------------------------------------------

def read_gencode_gtf(filepath):
    """Return (enst2ensg, enst2type) dicts from a Gencode GTF."""
    enst2ensg = {}
    enst2type = {}
    print(f"Reading in {filepath}", file=sys.stderr)
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            if len(tmp) < 9 or tmp[2] != "transcript":
                continue
            stuff = tmp[8]
            ensg_id = enst_id = gene_type = transcript_type = gene_name = None
            for s in stuff.split(";"):
                s = s.strip()
                m = re.search(r'gene_id "(.+?)"', s)
                if m:
                    ensg_id = m.group(1)
                m = re.search(r'transcript_id "(.+?)"', s)
                if m:
                    enst_id = m.group(1)
                m = re.search(r'transcript_type "(.+?)"', s)
                if m:
                    transcript_type = m.group(1)
                m = re.search(r'gene_type "(.+?)"', s)
                if m:
                    gene_type = m.group(1)
            if enst_id and ensg_id:
                enst2ensg[enst_id] = ensg_id
                if transcript_type:
                    enst2type[enst_id] = transcript_type
    return enst2ensg, enst2type


def read_gencode_tablebrowser(filepath, enst2ensg, enst2type):
    """
    Build gencode_features[chr][strand][bucket] = list of feature strings.
    Features encode exon/intron/UTR regions for each transcript.
    """
    gencode_features = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    print(f"reading in {filepath}", file=sys.stderr)
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            if len(tmp) < 11:
                continue
            enst = tmp[1]
            if enst == "name":
                continue
            chr_ = tmp[2]
            strand = tmp[3]
            txstart = int(tmp[4])
            txstop = int(tmp[5])
            cdsstart = int(tmp[6])
            cdsstop = int(tmp[7])
            starts = [int(x) for x in tmp[9].split(",") if x]
            stops = [int(x) for x in tmp[10].split(",") if x]

            transcript_type = enst2type.get(enst)
            if not transcript_type:
                print(f"error transcript_type {transcript_type} {enst}", file=sys.stderr)
                continue

            tmp_features = []

            if transcript_type == "protein_coding":
                for i in range(len(starts)):
                    s = starts[i]
                    e = stops[i]
                    if strand == "+":
                        if e < cdsstart:
                            tmp_features.append(f"{enst}|5utr|{s}-{e}")
                        elif s > cdsstop:
                            tmp_features.append(f"{enst}|3utr|{s}-{e}")
                        elif s > cdsstart and e < cdsstop:
                            tmp_features.append(f"{enst}|CDS|{s}-{e}")
                        else:
                            cds_s, cds_e = s, e
                            if s <= cdsstart <= e:
                                tmp_features.append(f"{enst}|5utr|{s}-{cdsstart}")
                                cds_s = cdsstart
                            if s <= cdsstop <= e:
                                tmp_features.append(f"{enst}|3utr|{cdsstop}-{e}")
                                cds_e = cdsstop
                            tmp_features.append(f"{enst}|CDS|{cds_s}-{cds_e}")
                    else:  # "-" strand
                        if e < cdsstart:
                            tmp_features.append(f"{enst}|3utr|{s}-{e}")
                        elif s > cdsstop:
                            tmp_features.append(f"{enst}|5utr|{s}-{e}")
                        elif s > cdsstart and e < cdsstop:
                            tmp_features.append(f"{enst}|CDS|{s}-{e}")
                        else:
                            cds_s, cds_e = s, e
                            if s <= cdsstart <= e:
                                tmp_features.append(f"{enst}|3utr|{s}-{cdsstart}")
                                cds_s = cdsstart
                            if s <= cdsstop <= e:
                                tmp_features.append(f"{enst}|5utr|{cdsstop}-{e}")
                                cds_e = cdsstop
                            tmp_features.append(f"{enst}|CDS|{cds_s}-{cds_e}")

                for i in range(len(starts) - 1):
                    gap = starts[i + 1] - stops[i]
                    if gap > 2 * 500:
                        tmp_features.append(f"{enst}|proxintron|{stops[i]}-{stops[i]+500}")
                        tmp_features.append(f"{enst}|distintron|{stops[i]+500}-{starts[i+1]-500}")
                        tmp_features.append(f"{enst}|proxintron|{starts[i+1]-500}-{starts[i+1]}")
                    else:
                        midpoint = (starts[i + 1] + stops[i]) // 2
                        tmp_features.append(f"{enst}|proxintron|{stops[i]}-{midpoint}")
                        tmp_features.append(f"{enst}|proxintron|{midpoint}-{starts[i+1]}")
            else:
                for i in range(len(starts)):
                    tmp_features.append(f"{enst}|noncoding_exon|{starts[i]}-{stops[i]}")
                for i in range(len(starts) - 1):
                    gap = starts[i + 1] - stops[i]
                    if gap > 2 * 500:
                        tmp_features.append(f"{enst}|noncoding_proxintron|{stops[i]}-{stops[i]+500}")
                        tmp_features.append(f"{enst}|noncoding_distintron|{stops[i]+500}-{starts[i+1]-500}")
                        tmp_features.append(f"{enst}|noncoding_proxintron|{starts[i+1]-500}-{starts[i+1]}")
                    else:
                        midpoint = (starts[i + 1] + stops[i]) // 2
                        tmp_features.append(f"{enst}|noncoding_proxintron|{stops[i]}-{midpoint}")
                        tmp_features.append(f"{enst}|noncoding_proxintron|{midpoint}-{starts[i+1]}")

            for feature in tmp_features:
                parts = feature.split("|")
                feat_enst, feat_type, feat_region = parts[0], parts[1], parts[2]
                reg_start, reg_stop = map(int, feat_region.split("-"))
                x = reg_start // GENOME_HASHING_VALUE
                y = reg_stop // GENOME_HASHING_VALUE
                for j in range(x, y + 1):
                    gencode_features[chr_][strand][j].append(feature)
                    gencode_features[chr_][REVSTRAND[strand]][j].append(
                        f"{feat_enst}|antisense_gencode|{reg_start}-{reg_stop}"
                    )

    return gencode_features


def read_peakfile(filepath):
    """
    Return peaks[chr][strand][bucket] = list of peak strings.
    BED format: chr start stop gene pval strand
    """
    peaks = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            if len(tmp) < 6:
                continue
            chr_ = tmp[0]
            start = int(tmp[1])
            stop = int(tmp[2])
            gene = tmp[3].upper()
            pval = tmp[4]
            strand = tmp[5]
            if chr_ == "genoName":
                continue
            x = start // GENOME_HASHING_VALUE
            y = stop // GENOME_HASHING_VALUE
            peak_sense = f"{chr_}:{start}-{stop}:{strand}:{gene}"
            peak_anti = f"{chr_}:{start}-{stop}:{strand}:antisense_{gene}"
            for i in range(x, y + 1):
                peaks[chr_][strand][i].append(peak_sense)
                peaks[chr_][REVSTRAND[strand]][i].append(peak_anti)
    return peaks


def read_filelists(filepath):
    """
    Return (enst2gene, convert_enst2type, convert_enst2priorityN)
    with uppercase ENST keys (matches Perl `uc($allenst)`).
    """
    enst2gene = {}
    convert_enst2type = {}
    convert_enst2priorityN = {}
    priority_n = 0
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            parts = line.split("\t")
            allenst = parts[0] if len(parts) > 0 else ""
            if not allenst:
                print(f"error missing enst {line} {filepath}", file=sys.stderr)
                continue
            allenst = allenst.upper()
            allensg = parts[1] if len(parts) > 1 else ""
            gid = parts[2] if len(parts) > 2 else ""
            type_label = parts[3] if len(parts) > 3 else ""
            gid = gid.rstrip("?").rstrip("_")
            type_label = type_label.rstrip("?")
            for enst in allenst.split("|"):
                enst2gene[enst] = gid
                enst2gene["antisense_" + enst] = "antisense_" + gid
                convert_enst2type[enst] = type_label
                convert_enst2type["antisense_" + enst] = "antisense_" + type_label
                convert_enst2priorityN[enst.lower()] = priority_n
                convert_enst2priorityN[(enst + "_DOUBLEMAP").lower()] = priority_n + 0.25
                convert_enst2priorityN[("antisense_" + enst).lower()] = priority_n + 0.5
                convert_enst2priorityN[("antisense_" + enst + "_DOUBLEMAP").lower()] = priority_n + 0.75
                priority_n += 1
    return enst2gene, convert_enst2type, convert_enst2priorityN


# ---------------------------------------------------------------------------
# CIGAR / alignment scoring
# ---------------------------------------------------------------------------

def parse_cigar_string(region_start_pos, flags, chr_, strand):
    """
    Convert SAM CIGAR to a list of interval strings:
        "{chr}:{strand}:{start}-{stop}"
    using 0-based half-open coordinates (matching Perl).
    """
    current_pos = region_start_pos
    region_start = region_start_pos
    regions = []
    for m in re.finditer(r"(\d+)([A-Z])", flags):
        length = int(m.group(1))
        op = m.group(2)
        if op == "N":
            regions.append(f"{chr_}:{strand}:{region_start - 1}-{current_pos - 1}")
            current_pos += length
            region_start = current_pos
        elif op == "M":
            current_pos += length
        elif op in ("S", "I"):
            pass  # soft-clip / insertion: no genome advance
        elif op == "D":
            current_pos += length
        # else: unknown op — ignore
    regions.append(f"{chr_}:{strand}:{region_start - 1}-{current_pos - 1}")
    return regions


def parse_cigar_for_alignment_score(cigar, phred):
    """
    Returns (gap_score, [insertion_strings]).
    insertion_strings are "{read_pos}|{length}".
    Simple case: pure match ("\\d+M") returns (0, []).
    """
    if re.match(r"^\d+M$", cigar):
        return 0, []

    gap_open = 5
    gap_extend = 3
    mm_penalty = 0
    insertions = []
    current_read_pos = 0

    for m in re.finditer(r"(\d+)([A-Z])", cigar):
        length = int(m.group(1))
        op = m.group(2)
        if op == "N":
            current_read_pos += length
        elif op == "M":
            current_read_pos += length
        elif op == "S":
            current_read_pos += length
            mm_penalty += gap_open + length * gap_extend
        elif op == "I":
            insertions.append(f"{current_read_pos}|{length}")
            current_read_pos += length
            mm_penalty += gap_open + length * gap_extend
        elif op == "D":
            mm_penalty += gap_open + length * gap_extend
    return mm_penalty, insertions


def parse_mismatch_string_for_score(md_string, phred, insertions, read_seq):
    """
    Compute mismatch penalty from MD string, quality scores, and insertions.
    Returns integer penalty.
    """
    if re.match(r"^\d+$", md_string):
        return 0

    quality_scores = [ord(c) - 33 for c in phred]
    mn, mx = 2, 6

    # Build insertion lookup: read_pos -> length
    ins_list = []
    for ins in insertions:
        pos_s, len_s = ins.split("|")
        ins_list.append((int(pos_s), int(len_s)))
    ins_idx = 0

    current_read_pos = 0
    mm_penalty = 0
    flags = md_string

    while flags:
        m = re.match(r"^(\d+)", flags)
        if m:
            advance = int(m.group(1))
            current_read_pos += advance
            # Skip over insertions that fall within this advance
            while ins_idx < len(ins_list) and current_read_pos >= ins_list[ins_idx][0]:
                current_read_pos += ins_list[ins_idx][1]
                ins_idx += 1
            flags = flags[len(m.group(1)):]
        elif re.match(r"^\^([A-Z]+)\d", flags):
            m2 = re.match(r"^\^([A-Z]+)", flags)
            length = len(m2.group(1))
            flags = flags[length + 1:]
        elif re.match(r"^([A-Z]+)", flags):
            m3 = re.match(r"^([A-Z]+)", flags)
            mismatch_bases = m3.group(1)
            for j in range(len(mismatch_bases)):
                pos = current_read_pos + j
                if pos >= len(quality_scores):
                    break
                base_score = quality_scores[pos]
                base = mismatch_bases[j]
                base_seq = read_seq[pos] if pos < len(read_seq) else "N"
                if base == "N" or base_seq == "N":
                    base_mm = 1
                elif (base in "ACGT") and (base_seq in "ACGT"):
                    base_mm = mn + math.floor((mx - mn) * (min(base_score, 40) / 40))
                else:
                    base_mm = 0  # unexpected
                mm_penalty += base_mm
            current_read_pos += len(mismatch_bases)
            flags = flags[len(mismatch_bases):]
        else:
            # Unrecognized — break to avoid infinite loop
            print(f"this is a flag I'm not expecting {flags}", file=sys.stderr)
            break

    return mm_penalty


def get_alignment_score(md_string, cigar, phred, read_seq):
    """Returns total alignment score (negative penalty, as in Perl)."""
    if md_string is None:
        return 0
    gap_score, insertions = parse_cigar_for_alignment_score(cigar, phred)
    mm_score = parse_mismatch_string_for_score(md_string, phred, insertions, read_seq)
    return 0 - mm_score - gap_score


def extract_md(flags_str):
    m = re.search(r"MD:Z:(\S+?)(\s|$)", flags_str)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# sort_priority helper (matches Perl `sort_priority`)
# ---------------------------------------------------------------------------

def sort_priority(enstlist_str, convert_enst2priorityN):
    ensts = enstlist_str.split("|")
    try:
        sorted_ensts = sorted(ensts, key=lambda e: convert_enst2priorityN.get(e.lower(), 0))
    except Exception:
        sorted_ensts = ensts
    return "|".join(sorted_ensts)


# ---------------------------------------------------------------------------
# Peak overlap helper
# ---------------------------------------------------------------------------

def find_overlapping_peaks(read_regions, peaks, chr_):
    """Return dict {peak_string: 1} for peaks overlapping any read region."""
    tmp_hash = {}
    for region in read_regions:
        rchr, rstr, rpos = region.split(":")
        rstart, rstop = map(int, rpos.split("-"))
        rx = rstart // GENOME_HASHING_VALUE
        ry = rstop // GENOME_HASHING_VALUE
        for ri in range(rx, ry + 1):
            for peak in peaks.get(rchr, {}).get(rstr, {}).get(ri, []):
                pchr, ppos, pstr, ptype = peak.split(":")
                pstart, pstop = map(int, ppos.split("-"))
                if pstart >= rstop or pstop <= rstart:
                    continue
                tmp_hash[peak] = 1
    return tmp_hash


def classify_genomic_read(
    chr_, frag_strand, read_start_position,
    peaks, chr_read, read_regions,
    gencode_features, enst2ensg, convert_enst2type, convert_enst2priorityN,
):
    """
    Determine ensttype and all_mapped_ensts for a uniquely-genomic read.
    Mirrors the SE and PE logic in read_unique_mapped_se/pe.
    """
    tmp_hash = find_overlapping_peaks(read_regions, peaks, chr_read)

    temp_peak_read_counts = {}
    converted_types = {}
    for peak in sorted(tmp_hash.keys()):  # sorted for determinism
        pchr, ppos, pstr, ptype = peak.split(":")
        temp_peak_read_counts[ptype + "_uniquegenomic"] = (
            temp_peak_read_counts.get(ptype + "_uniquegenomic", 0) + 1
        )
        if ptype not in convert_enst2type:
            print(f"peak {peak} {ptype}", file=sys.stderr)
        else:
            converted_types[convert_enst2type[ptype]] = 1

    sorted_types = sorted(temp_peak_read_counts.keys())
    all_mapped_ensts = "|".join(sorted_types)
    sorted_converted = sorted(converted_types.keys())
    all_ensttypes = "|".join(sorted_converted)

    if "miRNA" in converted_types or "miRNA-proximal" in converted_types:
        if "miRNA" in converted_types:
            final_feature_type = "miRNA"
        else:
            final_feature_type = "miRNA-proximal"
        all_ensttypes = "unique_" + final_feature_type
    elif len(temp_peak_read_counts) == 0:
        # No RepBase overlap — classify by Gencode feature
        rx = read_start_position // GENOME_HASHING_VALUE
        feature_flag = False
        tmp_gencode_hash = defaultdict(dict)
        for feat in gencode_features.get(chr_read, {}).get(frag_strand, {}).get(rx, []):
            gencode_enst, gencode_type, gencode_region = feat.split("|")
            gencode_start, gencode_stop = map(int, gencode_region.split("-"))
            if read_start_position < gencode_start:
                continue
            if read_start_position >= gencode_stop:
                continue
            gencode_ensg = enst2ensg.get(gencode_enst, gencode_enst)
            tmp_gencode_hash[gencode_type][gencode_ensg] = "contained"
            feature_flag = True

        final_feature_type = "intergenic"
        if feature_flag:
            if "CDS" in tmp_gencode_hash:
                final_feature_type = "CDS"
            elif "3utr" in tmp_gencode_hash or "5utr" in tmp_gencode_hash:
                if "3utr" in tmp_gencode_hash and "5utr" in tmp_gencode_hash:
                    final_feature_type = "5utr_and_3utr"
                elif "3utr" in tmp_gencode_hash:
                    final_feature_type = "3utr"
                else:
                    final_feature_type = "5utr"
            elif "proxintron" in tmp_gencode_hash:
                final_feature_type = "proxintron"
            elif "distintron" in tmp_gencode_hash:
                final_feature_type = "distintron"
            elif "noncoding_exon" in tmp_gencode_hash:
                final_feature_type = "noncoding_exon"
            elif "noncoding_proxintron" in tmp_gencode_hash:
                final_feature_type = "noncoding_proxintron"
            elif "noncoding_distintron" in tmp_gencode_hash:
                final_feature_type = "noncoding_distintron"
            elif "antisense_gencode" in tmp_gencode_hash:
                final_feature_type = "antisense_gencode"

        all_ensttypes = "unique_" + final_feature_type
        if final_feature_type in tmp_gencode_hash:
            sorted_keys = sorted(tmp_gencode_hash[final_feature_type].keys())
            all_mapped_ensts = "unique_" + "|".join(sorted_keys)
        else:
            all_mapped_ensts = "unique_" + final_feature_type

    return all_ensttypes, all_mapped_ensts


# ---------------------------------------------------------------------------
# Read rep-family SAM files
# ---------------------------------------------------------------------------

def read_rep_family_se(sam_file, read_hash, convert_enst2priorityN):
    with open(sam_file) as fh:
        for r1 in fh:
            r1 = r1.rstrip("\n")
            if r1.startswith("@"):
                continue
            tmp_r1 = r1.split("\t")
            r1name = tmp_r1[0].split()[0]
            r1sam_flag = int(tmp_r1[1])
            if r1sam_flag == 4:
                continue
            if r1sam_flag in (16, 272):
                frag_strand = "-"
            elif r1sam_flag in (0, 256):
                frag_strand = "+"
            else:
                continue

            flags_r1 = "\t".join(tmp_r1[11:])
            m = re.search(r"AS:i:(\S+?)(\s|$)", flags_r1)
            r1_score = int(m.group(1)) if m else 0
            total_score = r1_score

            m2 = re.search(r"ZZ:Z:(\S+?)(\s|$)", flags_r1)
            all_mapped_ensts = m2.group(1) if m2 else ""

            all_info = tmp_r1[2]
            all_ensttypes, all_primary_enst = all_info.split("||", 1) if "||" in all_info else (all_info, "")

            sorted_ensts = sort_priority(all_mapped_ensts, convert_enst2priorityN)
            read_hash[r1name] = {
                "rep_flag": 1,
                "rep_score": total_score,
                "file1flag": 1,
                "R1": r1,
                "mult_ensts": sorted_ensts,
                "ensttype": all_ensttypes,
            }


def read_rep_family_pe(sam_file, read_hash, convert_enst2priorityN):
    with open(sam_file) as fh:
        while True:
            r1 = fh.readline()
            if not r1:
                break
            r1 = r1.rstrip("\n")
            if r1.startswith("@"):
                continue
            r2 = fh.readline()
            if not r2:
                break
            r2 = r2.rstrip("\n")

            tmp_r1 = r1.split("\t")
            tmp_r2 = r2.split("\t")
            r1name = tmp_r1[0].split()[0]
            r2name = tmp_r2[0].split()[0]
            if r1name != r2name:
                print(f"paired end mismatch error: {sam_file} r1 {tmp_r1[0]} r2 {tmp_r2[0]}", file=sys.stderr)
            r1sam_flag = int(tmp_r1[1])
            if r1sam_flag == 77 or r1sam_flag == 141:
                continue

            if r1sam_flag in (99, 355):
                frag_strand = "-"
            elif r1sam_flag in (83, 339):
                frag_strand = "+"
            elif r1sam_flag in (147, 403):
                r1, r2 = r2, r1
                tmp_r1 = r1.split("\t")
                tmp_r2 = r2.split("\t")
                frag_strand = "-"
            elif r1sam_flag in (163, 419):
                r1, r2 = r2, r1
                tmp_r1 = r1.split("\t")
                tmp_r2 = r2.split("\t")
                frag_strand = "+"
            else:
                continue

            flags_r1 = "\t".join(tmp_r1[11:])
            flags_r2 = "\t".join(tmp_r2[11:])
            m1 = re.search(r"AS:i:(\S+?)(\s|$)", flags_r1)
            m2 = re.search(r"AS:i:(\S+?)(\s|$)", flags_r2)
            r1_score = int(m1.group(1)) if m1 else 0
            r2_score = int(m2.group(1)) if m2 else 0
            total_score = r1_score + r2_score

            mz = re.search(r"ZZ:Z:(\S+?)(\s|$)", flags_r1)
            all_mapped_ensts = mz.group(1) if mz else ""

            all_info = tmp_r1[2]
            all_ensttypes, all_primary_enst = all_info.split("||", 1) if "||" in all_info else (all_info, "")

            sorted_ensts = sort_priority(all_mapped_ensts, convert_enst2priorityN)
            read_hash[r1name] = {
                "rep_flag": 1,
                "rep_score": total_score,
                "file1flag": 1,
                "R1": r1,
                "R2": r2,
                "mult_ensts": sorted_ensts,
                "ensttype": all_ensttypes,
            }


# ---------------------------------------------------------------------------
# Read unique-genomic (rmrep) SAM files
# ---------------------------------------------------------------------------

def open_sam_stream(sam_file):
    """Open a .sam/.tmp directly, or a .bam through samtools.

    Returns (filehandle, process_or_None). The caller must pass the process to
    close_sam_stream() so a samtools failure is not read as an empty file --
    the same silent-failure that made a missing bowtie2 look like zero
    repeat-family reads (issue -475.16).
    """
    if sam_file.endswith(".sam") or sam_file.endswith(".tmp"):
        return open(sam_file), None
    if sam_file.endswith(".bam"):
        proc = subprocess.Popen(
            ["samtools", "view", "-h", sam_file],
            stdout=subprocess.PIPE,
            text=True,
        )
        return proc.stdout, proc
    raise ValueError(f"couldn't figure out format of {sam_file}")


def close_sam_stream(fh, proc, sam_file):
    """Close the stream and abort if samtools exited non-zero."""
    fh.close()
    if proc is not None:
        returncode = proc.wait()
        if returncode != 0:
            sys.exit(
                f"samtools view failed with exit code {returncode} on {sam_file}"
            )


def read_unique_mapped_se(sam_file, read_hash, peaks, gencode_features,
                          enst2ensg, convert_enst2type, convert_enst2priorityN):
    fi2_count = 0
    fh, proc = open_sam_stream(sam_file)

    for r1 in fh:
        r1 = r1.rstrip("\n")
        if r1.startswith("@"):
            continue
        fi2_count += 1
        if fi2_count % 100000 == 0:
            print(f"read {fi2_count}", file=sys.stderr)

        tmp_r1 = r1.split("\t")
        r1name = tmp_r1[0].split()[0]
        r1sam_flag = int(tmp_r1[1])
        if r1sam_flag == 4:
            continue

        if r1sam_flag in (16, 272):
            frag_strand = "-"
        elif r1sam_flag in (0, 256):
            frag_strand = "+"
        else:
            continue

        r1_cigar = tmp_r1[5]
        r1_chr = tmp_r1[2]
        r1_start = int(tmp_r1[3])
        flags_r1 = "\t".join(tmp_r1[11:])
        r1_mismatch = extract_md(flags_r1)
        r1_phred = tmp_r1[10]
        r1_seq = tmp_r1[9]
        r1_mmscore = get_alignment_score(r1_mismatch, r1_cigar, r1_phred, r1_seq)
        total_mmscore = r1_mmscore

        # Conflict resolution: unique genome beats rep only if score > rep_score + 24
        if r1name in read_hash and read_hash[r1name].get("file1flag") == 1:
            rep_score = read_hash[r1name].get("rep_score", 0)
            if total_mmscore > rep_score + 24:
                del read_hash[r1name]
            else:
                continue  # rep wins — discard unique genome mapping

        read_regions = parse_cigar_string(r1_start, r1_cigar, r1_chr, frag_strand)

        if frag_strand == "+":
            read_start_position = r1_start
        else:
            last_region = read_regions[-1]
            _, _, rpos = last_region.split(":")
            _, rstop = rpos.split("-")
            read_start_position = int(rstop) - 1

        all_ensttypes, all_mapped_ensts = classify_genomic_read(
            r1_chr, frag_strand, read_start_position,
            peaks, r1_chr, read_regions,
            gencode_features, enst2ensg, convert_enst2type, convert_enst2priorityN,
        )

        read_hash[r1name] = {
            "rep_flag": 0,
            "rep_score": total_mmscore,
            "file2flag": 1,
            "R1": r1,
            "mult_ensts": all_mapped_ensts,
            "ensttype": all_ensttypes,
        }

    close_sam_stream(fh, proc, sam_file)


def read_unique_mapped_pe(sam_file, read_hash, peaks, gencode_features,
                          enst2ensg, convert_enst2type, convert_enst2priorityN):
    fi2_count = 0
    fh, proc = open_sam_stream(sam_file)

    while True:
        r1_line = fh.readline()
        if not r1_line:
            break
        r1 = r1_line.rstrip("\n")
        if r1.startswith("@"):
            continue
        fi2_count += 1
        if fi2_count % 100000 == 0:
            print(f"read {fi2_count}", file=sys.stderr)
        r2 = fh.readline().rstrip("\n")

        tmp_r1 = r1.split("\t")
        tmp_r2 = r2.split("\t")
        r1name = tmp_r1[0].split()[0]
        r2name = tmp_r2[0].split()[0]
        if tmp_r1[0] != tmp_r2[0]:
            print(f"paired end mismatch error: {tmp_r1[0]} {tmp_r2[0]}", file=sys.stderr)

        r1sam_flag = int(tmp_r1[1])
        if r1sam_flag == 77 or r1sam_flag == 141:
            continue

        if r1sam_flag == 99:
            frag_strand = "-"
        elif r1sam_flag == 83:
            frag_strand = "+"
        elif r1sam_flag == 147:
            frag_strand = "-"
            r1, r2 = r2, r1
            tmp_r1, tmp_r2 = r1.split("\t"), r2.split("\t")
        elif r1sam_flag == 163:
            frag_strand = "+"
            r1, r2 = r2, r1
            tmp_r1, tmp_r2 = r1.split("\t"), r2.split("\t")
        else:
            continue

        r1_cigar = tmp_r1[5]
        r2_cigar = tmp_r2[5]
        r1_chr = tmp_r1[2]
        r1_start = int(tmp_r1[3])
        r2_chr = tmp_r2[2]
        r2_start = int(tmp_r2[3])

        flags_r1 = "\t".join(tmp_r1[11:])
        flags_r2 = "\t".join(tmp_r2[11:])
        r1_mismatch = extract_md(flags_r1)
        r2_mismatch = extract_md(flags_r2)
        r1_phred = tmp_r1[10]
        r2_phred = tmp_r2[10]
        r1_seq = tmp_r1[9]
        r2_seq = tmp_r2[9]
        r1_mmscore = get_alignment_score(r1_mismatch, r1_cigar, r1_phred, r1_seq)
        r2_mmscore = get_alignment_score(r2_mismatch, r2_cigar, r2_phred, r2_seq)
        total_mmscore = r1_mmscore + r2_mmscore

        if r1name in read_hash and read_hash[r1name].get("file1flag") == 1:
            rep_score = read_hash[r1name].get("rep_score", 0)
            if total_mmscore > rep_score + 24:
                del read_hash[r1name]
            else:
                continue

        # Use R2 for peak overlap (matches PE Perl code)
        read_regions = parse_cigar_string(r2_start, r2_cigar, r2_chr, frag_strand)

        if frag_strand == "+":
            read_start_position = r2_start
        else:
            last_region = read_regions[-1]
            _, _, rpos = last_region.split(":")
            _, rstop = rpos.split("-")
            read_start_position = int(rstop) - 1

        all_ensttypes, all_mapped_ensts = classify_genomic_read(
            r1_chr, frag_strand, read_start_position,
            peaks, r1_chr, read_regions,
            gencode_features, enst2ensg, convert_enst2type, convert_enst2priorityN,
        )

        read_hash[r1name] = {
            "rep_flag": 0,
            "rep_score": total_mmscore,
            "file2flag": 1,
            "R1": r1,
            "R2": r2,
            "mult_ensts": all_mapped_ensts,
            "ensttype": all_ensttypes,
        }

    close_sam_stream(fh, proc, sam_file)


# ---------------------------------------------------------------------------
# PCR duplicate removal
# ---------------------------------------------------------------------------

def run_pcr_duplicate_removal_se(read_hash, out_fh, predup_fh):
    """
    Returns (all_count, duplicate_count, unique_count_old,
             unique_genomic_count, unique_repfamily_count,
             count, count_enst, total_unique_mapped_read_num,
             rep_family_reads, unique_genomic).
    """
    fragment_hash = defaultdict(dict)
    all_count = duplicate_count = unique_count_old = 0
    unique_genomic_count = unique_repfamily_count = 0
    count = {}
    count_enst = {}
    total_unique_mapped_read_num = 0
    rep_family_reads = 0
    unique_genomic = 0

    for r1name in sorted(read_hash.keys()):  # sorted for determinism (matches Perl)
        r1 = read_hash[r1name]["R1"]
        tmp_r1 = r1.split("\t")
        r1sam_flag = int(tmp_r1[1])
        if r1sam_flag == 4:
            continue

        if r1sam_flag in (16, 272):
            frag_strand = "-"
        elif r1sam_flag in (0, 256):
            frag_strand = "+"
        else:
            continue

        read_name_parts = tmp_r1[0].split("_")
        randommer = read_name_parts[-1]

        r1_cigar = tmp_r1[5]
        r1_chr = tmp_r1[2]
        r1_start = int(tmp_r1[3])

        read_regions = parse_cigar_string(r1_start, r1_cigar, r1_chr, frag_strand)

        if frag_strand == "+":
            first_region = read_regions[-1]
            rchr, rstr, rpos = first_region.split(":")
            rstart, rstop = rpos.split("-")
            hashing_value = f"{rchr}:{rstr}:{rstart}\t{rchr}:{rstr}:{rstop}"
        else:
            first_region = read_regions[0]
            rchr, rstr, rpos = first_region.split(":")
            rstart, rstop = rpos.split("-")
            hashing_value = f"{rchr}:{rstr}:{rstart}\t{rchr}:{rstr}:{rstop}"

        hashing_value_toprint = f"{r1_chr}|{frag_strand}::{hashing_value}:{randommer}"
        predup_fh.write(f"{r1}\t{hashing_value_toprint}\n")

        all_count += 1
        frag_key = r1_chr + "|" + frag_strand
        umi_key = hashing_value + ":" + randommer

        if umi_key in fragment_hash[frag_key]:
            duplicate_count += 1
            del read_hash[r1name]
            continue
        else:
            fragment_hash[frag_key][umi_key] = 1

            if read_hash[r1name].get("file1flag") == 1:
                out_fh.write(f"{r1}\tRepFamily\t{r1_chr}\n")
                tmp_r1.append("RepFamily")
                tmp_r1.append(read_hash[r1name]["ensttype"] + "||" + read_hash[r1name]["mult_ensts"])
                unique_repfamily_count += 1
            elif read_hash[r1name].get("file2flag") == 1:
                out_fh.write(
                    f"{r1}\tUniqueGenomic\t"
                    f"{read_hash[r1name]['ensttype']}||{read_hash[r1name]['mult_ensts']}\n"
                )
                tmp_r1.append("UniqueGenomic")
                tmp_r1.append(read_hash[r1name]["ensttype"] + "||" + read_hash[r1name]["mult_ensts"])
                unique_genomic_count += 1
            else:
                print(
                    f"this shouldn't be hit - fi1flag {read_hash[r1name].get('file1flag')} "
                    f"fi2flag {read_hash[r1name].get('file2flag')} {r1name}",
                    file=sys.stderr,
                )

        unique_count_old += 1

        repmap_info = tmp_r1.pop()
        type_ = tmp_r1.pop()
        mult_transcripts = tmp_r1.pop()  # noqa: F841 (matches Perl pop but unused)

        ensttype, mult_ensts = repmap_info.split("||", 1) if "||" in repmap_info else (repmap_info, "")

        if type_ == "RepFamily":
            rep_family_reads += 1
        elif type_ == "UniqueGenomic":
            if r1_chr == "chrM":
                ensttype = f"chrM_unique_{frag_strand}strand"
            unique_genomic += 1

        if re.search(r"Simple_repeat", ensttype):
            ensttype = "Simple_repeat"

        count[ensttype] = count.get(ensttype, 0) + 1
        key = ensttype + "||" + mult_ensts
        count_enst[key] = count_enst.get(key, 0) + 1
        total_unique_mapped_read_num += 1

    return (all_count, duplicate_count, unique_count_old,
            unique_genomic_count, unique_repfamily_count,
            count, count_enst, total_unique_mapped_read_num,
            rep_family_reads, unique_genomic)


def run_pcr_duplicate_removal_pe(read_hash, out_fh, predup_fh):
    fragment_hash = defaultdict(dict)
    all_count = duplicate_count = unique_count_old = 0
    unique_genomic_count = unique_repfamily_count = 0
    count = {}
    count_enst = {}
    total_unique_mapped_read_num = 0
    rep_family_reads = 0
    unique_genomic = 0

    for r1name in sorted(read_hash.keys()):  # sorted for determinism
        r1 = read_hash[r1name]["R1"]
        r2 = read_hash[r1name]["R2"]
        tmp_r1 = r1.split("\t")
        tmp_r2 = r2.split("\t")
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

        randommer = tmp_r1[0].split(":")[0]

        r1_cigar = tmp_r1[5]
        r2_cigar = tmp_r2[5]
        r1_chr = tmp_r1[2]
        r1_start = int(tmp_r1[3])
        r2_chr = tmp_r2[2]
        r2_start = int(tmp_r2[3])

        read_regions_r1 = parse_cigar_string(r1_start, r1_cigar, r1_chr, frag_strand)
        read_regions_r2 = parse_cigar_string(r2_start, r2_cigar, r2_chr, frag_strand)

        if frag_strand == "+":
            r1_fr = read_regions_r1[-1]
            r1_rchr, r1_rstr, r1_rpos = r1_fr.split(":")
            r1_rstart, r1_rstop = r1_rpos.split("-")
            r2_fr = read_regions_r2[0]
            r2_rchr, r2_rstr, r2_rpos = r2_fr.split(":")
            r2_rstart, r2_rstop = r2_rpos.split("-")
            hashing_value = (
                f"{r2_rchr}:{r2_rstr}:{r2_rstart}\t{r1_rchr}:{r1_rstr}:{r1_rstop}"
            )
        else:
            r1_fr = read_regions_r1[0]
            r1_rchr, r1_rstr, r1_rpos = r1_fr.split(":")
            r1_rstart, r1_rstop = r1_rpos.split("-")
            r2_fr = read_regions_r2[-1]
            r2_rchr, r2_rstr, r2_rpos = r2_fr.split(":")
            r2_rstart, r2_rstop = r2_rpos.split("-")
            hashing_value = (
                f"{r1_rchr}:{r1_rstr}:{r1_rstart}\t{r2_rchr}:{r2_rstr}:{r2_rstop}"
            )

        hashing_value_toprint = f"{r1_chr}|{frag_strand}::{hashing_value}:{randommer}"
        predup_fh.write(f"{r1}\t{hashing_value_toprint}\n{r2}\t{hashing_value_toprint}\n")

        all_count += 1
        frag_key = r1_chr + "|" + frag_strand
        umi_key = hashing_value + ":" + randommer

        if umi_key in fragment_hash[frag_key]:
            duplicate_count += 1
            del read_hash[r1name]
            continue
        else:
            fragment_hash[frag_key][umi_key] = 1

            if read_hash[r1name].get("file1flag") == 1:
                out_fh.write(f"{r1}\tRepFamily\t{r1_chr}\n{r2}\tRepFamily\t{r1_chr}\n")
                tmp_r1.append("RepFamily")
                tmp_r1.append(read_hash[r1name]["ensttype"] + "||" + read_hash[r1name]["mult_ensts"])
                unique_repfamily_count += 1
            elif read_hash[r1name].get("file2flag") == 1:
                out_fh.write(
                    f"{r1}\tUniqueGenomic\t"
                    f"{read_hash[r1name]['ensttype']}||{read_hash[r1name]['mult_ensts']}\n"
                    f"{r2}\tUniqueGenomic\t"
                    f"{read_hash[r1name]['ensttype']}||{read_hash[r1name]['mult_ensts']}\n"
                )
                tmp_r1.append("UniqueGenomic")
                tmp_r1.append(read_hash[r1name]["ensttype"] + "||" + read_hash[r1name]["mult_ensts"])
                unique_genomic_count += 1
            else:
                print(
                    f"this shouldn't be hit - fi1flag {read_hash[r1name].get('file1flag')} "
                    f"fi2flag {read_hash[r1name].get('file2flag')} {r1name}",
                    file=sys.stderr,
                )

        unique_count_old += 1

        repmap_info = tmp_r1.pop()
        type_ = tmp_r1.pop()
        mult_transcripts = tmp_r1.pop()  # noqa: F841

        ensttype, mult_ensts = repmap_info.split("||", 1) if "||" in repmap_info else (repmap_info, "")

        if type_ == "RepFamily":
            rep_family_reads += 1
        elif type_ == "UniqueGenomic":
            if r1_chr == "chrM":
                ensttype = f"chrM_unique_{frag_strand}strand"
            unique_genomic += 1

        if re.search(r"Simple_repeat", ensttype):
            ensttype = "Simple_repeat"

        count[ensttype] = count.get(ensttype, 0) + 1
        key = ensttype + "||" + mult_ensts
        count_enst[key] = count_enst.get(key, 0) + 1
        total_unique_mapped_read_num += 1

    return (all_count, duplicate_count, unique_count_old,
            unique_genomic_count, unique_repfamily_count,
            count, count_enst, total_unique_mapped_read_num,
            rep_family_reads, unique_genomic)


# ---------------------------------------------------------------------------
# Output: write parsed_v2 file
# ---------------------------------------------------------------------------

def write_parsed(
    count_out, enst2gene,
    all_count, duplicate_count, unique_count_old,
    unique_genomic_count, unique_repfamily_count,
    count, count_enst, total_unique_mapped_read_num,
    rep_family_reads, unique_genomic,
):
    with open(count_out, "w") as cnt:
        cnt.write(
            f"#READINFO\tAll reads:\t{all_count}\tPCR duplicates removed:\t{duplicate_count}"
            f"\tUsable Remaining:\t{unique_count_old}\tUsable from genomic mapping:\t{unique_genomic_count}"
            f"\tUsable from family mapping:\t{unique_repfamily_count}\n"
        )
        cnt.write(f"#READINFO\tUsableReads\t{total_unique_mapped_read_num}\n")
        if total_unique_mapped_read_num > 0:
            cnt.write(
                f"#READINFO\tGenomicReads\t{unique_genomic}\t"
                f"{unique_genomic / total_unique_mapped_read_num:.5f}\n"
            )
            cnt.write(
                f"#READINFO\tRepFamilyReads\t{rep_family_reads}\t"
                f"{rep_family_reads / total_unique_mapped_read_num:.5f}\n"
            )
        else:
            cnt.write(f"#READINFO\tGenomicReads\t{unique_genomic}\t0\n")
            cnt.write(f"#READINFO\tRepFamilyReads\t{rep_family_reads}\t0\n")

        # TOTAL lines sorted by count descending
        for k in sorted(count.keys(), key=lambda x: -count[x]):
            rpm = count[k] * 1000000 / total_unique_mapped_read_num if total_unique_mapped_read_num > 0 else 0
            cnt.write(f"TOTAL\t{k}\t{count[k]}\t{rpm:.5f}\n")

        # ELEMENT lines sorted by count descending
        for s in sorted(count_enst.keys(), key=lambda x: -count_enst[x]):
            ensttype, multensts = s.split("||", 1) if "||" in s else (s, "")
            multensts_groups = multensts.split("|")
            genes_final = []
            type_ = ensttype
            for group in multensts_groups:
                gids = group.split(";;")
                genes_short = []
                for gid in gids:
                    genes_short.append(enst2gene.get(gid, gid))
                genes_final.append(";;".join(genes_short))
            ensg_all = "|".join(genes_final)
            readnum = count_enst[s]
            rpm = readnum * 1000000 / total_unique_mapped_read_num if total_unique_mapped_read_num > 0 else 0
            cnt.write(
                f"ELEMENT\t{type_}\t{readnum}\t{rpm:.5f}\t{s}\t{ensg_all}\n"
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 8:
        print(
            "Usage: deduplicate.py <rep_tmp> <rmrep_tmp> <SE_or_PE> "
            "<gencode_gtf> <gencode_table_browser> <rep_mask_bed> <file_list>",
            file=sys.stderr,
        )
        sys.exit(1)

    repfamily_sam = sys.argv[1]
    gabe_rmrep_sam = sys.argv[2]
    se_or_pe = sys.argv[3]
    gencode_gtf_file = sys.argv[4]
    gencode_tablebrowser_file = sys.argv[5]
    repmask_bed = sys.argv[6]
    filelist_file = sys.argv[7]

    repfamily_sam_short = os.path.basename(repfamily_sam)

    # Load reference data
    enst2ensg, enst2type = read_gencode_gtf(gencode_gtf_file)
    gencode_features = read_gencode_tablebrowser(gencode_tablebrowser_file, enst2ensg, enst2type)
    peaks = read_peakfile(repmask_bed)
    enst2gene, convert_enst2type, convert_enst2priorityN = read_filelists(filelist_file)

    # Output file names (CWD, from basename of rep_tmp)
    output_fi = repfamily_sam_short + ".combined_w_uniquemap.rmDup.sam"
    pre_rmdup_fi = repfamily_sam_short + ".combined_w_uniquemap.prermDup.sam"
    count_out = output_fi + ".parsed_v2.20201210.txt"
    out_done = count_out + ".done"

    read_hash = {}

    with open(output_fi, "w") as out_fh, open(pre_rmdup_fi, "w") as predup_fh:

        if se_or_pe == "SE":
            print("running in single-end mode", file=sys.stderr)
            read_rep_family_se(repfamily_sam, read_hash, convert_enst2priorityN)
            read_unique_mapped_se(
                gabe_rmrep_sam, read_hash, peaks, gencode_features,
                enst2ensg, convert_enst2type, convert_enst2priorityN,
            )
            (all_count, duplicate_count, unique_count_old,
             unique_genomic_count, unique_repfamily_count,
             count, count_enst, total_unique_mapped_read_num,
             rep_family_reads, unique_genomic) = run_pcr_duplicate_removal_se(
                read_hash, out_fh, predup_fh
            )

        elif se_or_pe == "PE":
            print("running in paired-end mode", file=sys.stderr)
            read_rep_family_pe(repfamily_sam, read_hash, convert_enst2priorityN)
            read_unique_mapped_pe(
                gabe_rmrep_sam, read_hash, peaks, gencode_features,
                enst2ensg, convert_enst2type, convert_enst2priorityN,
            )
            (all_count, duplicate_count, unique_count_old,
             unique_genomic_count, unique_repfamily_count,
             count, count_enst, total_unique_mapped_read_num,
             rep_family_reads, unique_genomic) = run_pcr_duplicate_removal_pe(
                read_hash, out_fh, predup_fh
            )
        else:
            print(f"Unknown SE_or_PE value: {se_or_pe}", file=sys.stderr)
            sys.exit(1)

    write_parsed(
        count_out, enst2gene,
        all_count, duplicate_count, unique_count_old,
        unique_genomic_count, unique_repfamily_count,
        count, count_enst, total_unique_mapped_read_num,
        rep_family_reads, unique_genomic,
    )

    with open(out_done, "w") as donef:
        donef.write("jobs done\n")


if __name__ == "__main__":
    main()
