#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass
class ReadEntry:
    read_name: str
    r1: str
    r2: Optional[str]
    source_type: str
    ensttype: str
    mult_ensts: str
    strand: str
    chrom: str
    umi: str
    hash_key: str


def open_maybe_gzip(path: Path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return path.open('r', encoding='utf-8', errors='replace')


def parse_cigar_regions(start: int, cigar: str, chrom: str, strand: str) -> List[Tuple[int, int]]:
    pos = start
    regions: List[Tuple[int, int]] = []
    region_start = pos
    for length_s, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length_s)
        if op in {'M', '=', 'X'}:
            pos += length
        elif op in {'D', 'N'}:
            if region_start < pos:
                regions.append((region_start, pos - 1))
            pos += length
            region_start = pos
        else:
            # I/S/H/P do not consume reference.
            pass
    if region_start < pos:
        regions.append((region_start, pos - 1))
    if not regions:
        regions.append((start, start))
    return regions


def choose_strand_se(flag: int) -> Optional[str]:
    if flag in (16, 272):
        return '-'
    if flag in (0, 256):
        return '+'
    return None


def choose_strand_pe(flag: int) -> Optional[str]:
    if flag in (99, 355, 147, 403):
        return '-'
    if flag in (83, 339, 163, 419):
        return '+'
    return None


def parse_umi_se(qname: str) -> str:
    parts = qname.split('_')
    return parts[-1] if parts else qname


def parse_umi_pe(qname: str) -> str:
    parts = qname.split(':')
    return parts[0] if parts else qname


def make_hash_key_se(cols: List[str], strand: str, umi: str) -> str:
    chrom = cols[2]
    start = int(cols[3])
    regions = parse_cigar_regions(start, cols[5], chrom, strand)
    if strand == '+':
        first = regions[-1]
    else:
        first = regions[0]
    hashing_value = f'{chrom}:{strand}:{first[0]}\t{chrom}:{strand}:{first[1]}'
    return f'{hashing_value}:{umi}'


def make_hash_key_pe(cols1: List[str], cols2: List[str], strand: str, umi: str) -> str:
    c1 = cols1[2]
    s1 = int(cols1[3])
    r1 = parse_cigar_regions(s1, cols1[5], c1, strand)
    c2 = cols2[2]
    s2 = int(cols2[3])
    r2 = parse_cigar_regions(s2, cols2[5], c2, strand)
    if strand == '+':
        a = r1[-1]
        b = r2[0]
    else:
        a = r1[0]
        b = r2[-1]
    return f'{c1}:{strand}:{a[0]}\t{c1}:{strand}:{a[1]}::{c2}:{strand}:{b[0]}\t{c2}:{strand}:{b[1]}:{umi}'


def parse_rep_family_entries(path: Path, se_or_pe: str) -> Dict[str, ReadEntry]:
    entries: Dict[str, ReadEntry] = {}
    with open_maybe_gzip(path) as fh:
        if se_or_pe == 'SE':
            for raw in fh:
                if raw.startswith('@'):
                    continue
                line = raw.rstrip('\n')
                cols = line.split('\t')
                if len(cols) < 6:
                    continue
                flag = int(cols[1])
                if flag == 4:
                    continue
                strand = choose_strand_se(flag)
                if strand is None:
                    continue
                qname = cols[0]
                read_name = '_'.join(qname.split('_')[:-1]) if '_' in qname else qname
                umi = parse_umi_se(qname)
                repinfo = cols[2]
                ensttype = repinfo.split('||', 1)[0]
                zz = ''
                for tag in cols[11:]:
                    if tag.startswith('ZZ:Z:'):
                        zz = tag[5:]
                        break
                hash_key = make_hash_key_se(cols, strand, umi)
                entries[read_name] = ReadEntry(
                    read_name=read_name,
                    r1=line,
                    r2=None,
                    source_type='RepFamily',
                    ensttype=ensttype,
                    mult_ensts=zz if zz else repinfo,
                    strand=strand,
                    chrom=cols[2],
                    umi=umi,
                    hash_key=hash_key,
                )
        else:
            it = iter(fh)
            for raw1 in it:
                if raw1.startswith('@'):
                    continue
                try:
                    raw2 = next(it)
                except StopIteration:
                    break
                line1 = raw1.rstrip('\n')
                line2 = raw2.rstrip('\n')
                cols1 = line1.split('\t')
                cols2 = line2.split('\t')
                if len(cols1) < 6 or len(cols2) < 6:
                    continue
                flag = int(cols1[1])
                if flag in (77, 141):
                    continue
                strand = choose_strand_pe(flag)
                if strand is None:
                    continue
                read_name = cols1[0].split()[0]
                if read_name != cols2[0].split()[0]:
                    continue
                umi = parse_umi_pe(cols1[0])
                repinfo = cols1[2]
                ensttype = repinfo.split('||', 1)[0]
                zz = ''
                for tag in cols1[11:]:
                    if tag.startswith('ZZ:Z:'):
                        zz = tag[5:]
                        break
                hash_key = make_hash_key_pe(cols1, cols2, strand, umi)
                entries[read_name] = ReadEntry(
                    read_name=read_name,
                    r1=line1,
                    r2=line2,
                    source_type='RepFamily',
                    ensttype=ensttype,
                    mult_ensts=zz if zz else repinfo,
                    strand=strand,
                    chrom=cols1[2],
                    umi=umi,
                    hash_key=hash_key,
                )
    return entries


def parse_unique_genomic_entries(path: Path, se_or_pe: str) -> Dict[str, ReadEntry]:
    entries: Dict[str, ReadEntry] = {}
    with open_maybe_gzip(path) as fh:
        if se_or_pe == 'SE':
            for raw in fh:
                if raw.startswith('@'):
                    continue
                line = raw.rstrip('\n')
                cols = line.split('\t')
                if len(cols) < 6:
                    continue
                flag = int(cols[1])
                if flag == 4:
                    continue
                strand = choose_strand_se(flag)
                if strand is None:
                    continue
                qname = cols[0]
                read_name = '_'.join(qname.split('_')[:-1]) if '_' in qname else qname
                umi = parse_umi_se(qname)
                hash_key = make_hash_key_se(cols, strand, umi)
                entries[read_name] = ReadEntry(
                    read_name=read_name,
                    r1=line,
                    r2=None,
                    source_type='UniqueGenomic',
                    ensttype=cols[2],
                    mult_ensts=cols[2],
                    strand=strand,
                    chrom=cols[2],
                    umi=umi,
                    hash_key=hash_key,
                )
        else:
            it = iter(fh)
            for raw1 in it:
                if raw1.startswith('@'):
                    continue
                try:
                    raw2 = next(it)
                except StopIteration:
                    break
                line1 = raw1.rstrip('\n')
                line2 = raw2.rstrip('\n')
                cols1 = line1.split('\t')
                cols2 = line2.split('\t')
                if len(cols1) < 6 or len(cols2) < 6:
                    continue
                flag = int(cols1[1])
                if flag in (77, 141):
                    continue
                strand = choose_strand_pe(flag)
                if strand is None:
                    continue
                read_name = cols1[0].split()[0]
                if read_name != cols2[0].split()[0]:
                    continue
                umi = parse_umi_pe(cols1[0])
                hash_key = make_hash_key_pe(cols1, cols2, strand, umi)
                entries[read_name] = ReadEntry(
                    read_name=read_name,
                    r1=line1,
                    r2=line2,
                    source_type='UniqueGenomic',
                    ensttype=cols1[2],
                    mult_ensts=cols1[2],
                    strand=strand,
                    chrom=cols1[2],
                    umi=umi,
                    hash_key=hash_key,
                )
    return entries


def sorted_entries(entries: Dict[str, ReadEntry]) -> Iterable[ReadEntry]:
    for k in sorted(entries.keys()):
        yield entries[k]


def normalize_ensttype(ensttype: str, source_type: str, chrom: str, strand: str) -> str:
    out = ensttype
    if source_type == 'UniqueGenomic' and chrom == 'chrM':
        out = f'chrM_unique_{strand}strand'
    if 'Simple_repeat' in out:
        out = out.replace('Simple_repeat', 'Simple_repeat')
        out = out.replace('Simple_repeat', 'Simple_repeat')
    if 'Simple_repeat' not in out and 'Simple_repeat' in out.replace('Simple_repeat', 'Simple_repeat'):
        out = out.replace('Simple_repeat', 'Simple_repeat')
    if 'Simple_repeat' not in out and 'Simple_repeat' in ensttype:
        out = 'Simple_repeat'
    return out


def write_outputs(
    repfamily_sam: Path,
    combined_entries: Dict[str, ReadEntry],
) -> None:
    rep_short = repfamily_sam.name
    out_rmdup = Path(rep_short + '.combined_w_uniquemap.rmDup.sam')
    out_prermdup = Path(rep_short + '.combined_w_uniquemap.prermDup.sam')
    out_parsed = Path(str(out_rmdup) + '.parsed_v2.20201210.txt')
    out_done = Path(str(out_parsed) + '.done')

    all_count = 0
    duplicate_count = 0
    unique_count = 0
    unique_genomic_count = 0
    unique_repfamily_count = 0

    total_unique_mapped = 0
    rep_family_reads = 0
    unique_genomic_reads = 0

    total_count: Dict[str, int] = {}
    element_count: Dict[Tuple[str, str], int] = {}

    fragment_hash: Dict[Tuple[str, str], set[str]] = {}

    with out_prermdup.open('w', encoding='utf-8') as predup, out_rmdup.open('w', encoding='utf-8') as out:
        for entry in sorted_entries(combined_entries):
            all_count += 1
            chrom_key = entry.chrom
            strand_key = entry.strand
            full_key = entry.hash_key
            key_group = (chrom_key, strand_key)
            if key_group not in fragment_hash:
                fragment_hash[key_group] = set()

            predup_line = entry.r1 if entry.r2 is None else f'{entry.r1}\n{entry.r2}'
            predup.write(predup_line + '\t' + f'{entry.chrom}|{entry.strand}::{entry.hash_key}' + '\n')

            if full_key in fragment_hash[key_group]:
                duplicate_count += 1
                continue
            fragment_hash[key_group].add(full_key)

            unique_count += 1
            ensttype = normalize_ensttype(entry.ensttype, entry.source_type, entry.chrom, entry.strand)
            repmap_info = f'{ensttype}||{entry.mult_ensts}'
            if entry.source_type == 'RepFamily':
                unique_repfamily_count += 1
                rep_family_reads += 1
            else:
                unique_genomic_count += 1
                unique_genomic_reads += 1

            if entry.r2 is None:
                out.write(f'{entry.r1}\t{entry.source_type}\t{repmap_info}\n')
            else:
                out.write(f'{entry.r1}\t{entry.source_type}\t{repmap_info}\n')
                out.write(f'{entry.r2}\t{entry.source_type}\t{repmap_info}\n')

            total_unique_mapped += 1
            total_count[ensttype] = total_count.get(ensttype, 0) + 1
            ek = (f'{ensttype}||{entry.mult_ensts}', ensttype)
            element_count[ek] = element_count.get(ek, 0) + 1

    with out_parsed.open('w', encoding='utf-8') as count:
        count.write(
            '#READINFO\tAll reads:\t'
            + str(all_count)
            + '\tPCR duplicates removed:\t'
            + str(duplicate_count)
            + '\tUsable Remaining:\t'
            + str(unique_count)
            + '\tUsable from genomic mapping:\t'
            + str(unique_genomic_count)
            + '\tUsable from family mapping:\t'
            + str(unique_repfamily_count)
            + '\n'
        )
        count.write(f'#READINFO\tUsableReads\t{total_unique_mapped}\n')
        if total_unique_mapped > 0:
            count.write(f'#READINFO\tGenomicReads\t{unique_genomic_reads}\t{unique_genomic_reads/total_unique_mapped:.5f}\n')
            count.write(f'#READINFO\tRepFamilyReads\t{rep_family_reads}\t{rep_family_reads/total_unique_mapped:.5f}\n')
        else:
            count.write('#READINFO\tGenomicReads\t0\t0\n')
            count.write('#READINFO\tRepFamilyReads\t0\t0\n')

        for element, readnum in sorted(total_count.items(), key=lambda x: x[1], reverse=True):
            rpm = (readnum * 1000000 / total_unique_mapped) if total_unique_mapped else 0
            count.write(f'TOTAL\t{element}\t{readnum}\t{rpm:.5f}\n')

        for (enst_all, ensg_primary), readnum in sorted(element_count.items(), key=lambda x: x[1], reverse=True):
            rpm = (readnum * 1000000 / total_unique_mapped) if total_unique_mapped else 0
            # Keep column count/shape aligned with legacy parser consumers.
            count.write(f'ELEMENT\t{ensg_primary}\t{readnum}\t{rpm:.5f}\t{enst_all}\t{ensg_primary}\n')

    out_done.write_text('jobs done\n', encoding='utf-8')


def main() -> int:
    ap = argparse.ArgumentParser(description='Python replacement for duplicate_removal_inline_paired...perl')
    ap.add_argument('repFamilySam')
    ap.add_argument('rmRepSam')
    ap.add_argument('se_or_pe')
    ap.add_argument('gencodeGTF')
    ap.add_argument('gencodeTableBrowser')
    ap.add_argument('repMaskBedFile')
    ap.add_argument('fileList1')
    args = ap.parse_args()

    se_or_pe = args.se_or_pe.upper()
    if se_or_pe not in {'SE', 'PE'}:
        print('fatal error - SE or PE not defined', file=sys.stderr)
        return 2

    repfamily_sam = Path(args.repFamilySam)
    rmrep_sam = Path(args.rmRepSam)

    rep_entries = parse_rep_family_entries(repfamily_sam, se_or_pe)
    uniq_entries = parse_unique_genomic_entries(rmrep_sam, se_or_pe)

    combined = dict(rep_entries)
    for read_name, entry in uniq_entries.items():
        if read_name not in combined:
            combined[read_name] = entry

    write_outputs(repfamily_sam, combined)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
