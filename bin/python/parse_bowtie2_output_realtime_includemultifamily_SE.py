#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


RRNA_EXTRA_HASH = {
    'RNA28S': 'RNA45S',
    'RNA18S': 'RNA45S',
    'RNA5-8S': 'RNA45S',
    'antisense_RNA28S': 'antisense_RNA45S',
    'antisense_RNA18S': 'antisense_RNA45S',
    'antisense_RNA5-8S': 'antisense_RNA45S',
}
RRNA_EXTRA_HASH_REV = {
    'RNA45S': {'RNA28S', 'RNA18S', 'RNA5-8S'},
    'antisense_RNA45S': {'antisense_RNA28S', 'antisense_RNA18S', 'antisense_RNA5-8S'},
}


@dataclass
class ReadRecord:
    r1: Dict[str, str] = field(default_factory=dict)
    flags: Dict[str, int] = field(default_factory=dict)
    quality: int = 0
    mult_ensts: Dict[str, List[str]] = field(default_factory=lambda: defaultdict(list))
    enst: Dict[str, str] = field(default_factory=dict)
    master_enst: Dict[str, str] = field(default_factory=dict)


def read_in_filelists(filelist_file: Path) -> tuple[dict[str, str], dict[str, str]]:
    enst2gene: dict[str, str] = {}
    convert_enst2type: dict[str, str] = {}
    priority_n = 0
    with filelist_file.open('r', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            allenst, _allensg, gid, type_label = parts[:4]
            type_label = type_label.rstrip('_')
            if not allenst:
                print(f'error missing enst {line} {filelist_file}', file=sys.stderr)
                continue
            for enst in allenst.split('|'):
                enst2gene[enst] = gid
                convert_enst2type[enst] = f'{type_label}:{priority_n}'
                priority_n += 1
    return enst2gene, convert_enst2type


def parse_as_score(tags: list[str]) -> int:
    for tag in tags:
        if tag.startswith('AS:i:'):
            try:
                return int(tag[5:])
            except ValueError:
                pass
    return -10**9


def print_output(read_hash: dict[str, ReadRecord], samout, multimapping_hash: dict[str, int]) -> None:
    for read_name, rec in read_hash.items():
        ensttype_array = sorted(rec.flags.keys())
        if not ensttype_array:
            continue
        ensttype = ensttype_array[0]
        masterenst_array = [rec.master_enst[t] for t in ensttype_array if t in rec.master_enst]
        ensttype_join = '|'.join(ensttype_array)
        masterenst_join = '|'.join(masterenst_array)

        if len(rec.flags) == 1:
            r1_cols = rec.r1[ensttype].split('\t')
            r1_cols[2] = f'{ensttype_join}||{masterenst_join}'
            r1_line = '\t'.join(r1_cols)
            zz = '|'.join(rec.mult_ensts.get(ensttype, []))
            samout.write(f'{r1_line}\tZZ:Z:{zz}\n')
        else:
            all_mult_ensts: list[str] = []
            for key in ensttype_array:
                all_mult_ensts.append('|'.join(rec.mult_ensts.get(key, [])))
            final_mult_ensts = '|'.join(all_mult_ensts)
            if ensttype not in rec.r1:
                print(f"weird error - {read_name} {ensttype} readhash doesn't exist ? {rec.flags.get(ensttype)}", file=sys.stderr)
                continue
            r1_cols = rec.r1[ensttype].split('\t')
            r1_cols[2] = f'{ensttype_join}||{masterenst_join}'
            r1_line = '\t'.join(r1_cols)
            samout.write(f'{r1_line}\tZZ:Z:{final_mult_ensts}\n')
            multimapping_type = '|'.join(rec.flags.keys())
            multimapping_hash[multimapping_type] = multimapping_hash.get(multimapping_type, 0) + 1


def stream_sam_lines(args: argparse.Namespace):
    if args.sam_input:
        with Path(args.sam_input).open('r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                yield line
        return

    bowtie_out = f'{args.output}.bowtieout'
    cmd = [
        'bowtie2',
        '-q',
        '--sensitive',
        '-a',
        '-p',
        str(args.threads),
        '--no-mixed',
        '--reorder',
        '-x',
        args.bowtie_db,
        '-U',
        args.fastq_file1,
    ]
    print('command ' + ' '.join(cmd) + f' 2> {bowtie_out}', file=sys.stderr)
    with Path(bowtie_out).open('w', encoding='utf-8') as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(f'bowtie2 failed with exit code {ret}')


def main() -> int:
    ap = argparse.ArgumentParser(description='Python port of parse_bowtie2_output_realtime_includemultifamily_SE.pl')
    ap.add_argument('fastq_file1')
    ap.add_argument('bowtie_db')
    ap.add_argument('output')
    ap.add_argument('filelist_file')
    ap.add_argument('--sam-input', default=None, help='Optional SAM input stream file for testing/debugging')
    ap.add_argument('--threads', type=int, default=3)
    args = ap.parse_args()

    _enst2gene, convert_enst2type = read_in_filelists(Path(args.filelist_file))

    output_path = Path(args.output)
    multimapping_out = Path(args.output + '.multimapping_deleted')
    done_file = Path(args.output + '.done')

    read_hash: dict[str, ReadRecord] = {}
    multimapping_hash: dict[str, int] = {}
    print_batch = 10000
    read_counter = 0
    prev_r1name = ''

    with output_path.open('w', encoding='utf-8') as samout:
        for raw in stream_sam_lines(args):
            r1 = raw.rstrip('\n')
            if r1.startswith('@'):
                samout.write(r1 + '\n')
                continue

            cols = r1.split('\t')
            if len(cols) < 12:
                continue

            read_name_parts = cols[0].split('_')
            if len(read_name_parts) > 1:
                r1name = '_'.join(read_name_parts[:-1])
            else:
                r1name = cols[0]

            try:
                r1sam_flag = int(cols[1])
            except ValueError:
                continue
            if r1sam_flag == 4:
                continue

            if r1sam_flag in (16, 272):
                frag_strand = '-'
            elif r1sam_flag in (0, 256):
                frag_strand = '+'
            else:
                continue

            paired_mismatch_score = parse_as_score(cols[11:])

            mapped_enst = cols[2]
            mapped_enst_full = cols[2]
            if mapped_enst.endswith('_spliced'):
                mapped_enst = mapped_enst[: -len('_spliced')]
            if mapped_enst.endswith('_withgenomeflank'):
                mapped_enst = mapped_enst[: -len('_withgenomeflank')]

            if mapped_enst not in convert_enst2type:
                print(f'enst2type is missing for {mapped_enst} {r1}', file=sys.stderr)
                continue
            ensttype, enstpriority_s = convert_enst2type[mapped_enst].split(':', 1)
            enstpriority = int(enstpriority_s)

            if frag_strand == '-':
                ensttype = 'antisense_' + ensttype
                mapped_enst_full = 'antisense_' + mapped_enst_full

            if r1name != prev_r1name:
                read_counter += 1
                if read_counter > print_batch:
                    print_output(read_hash, samout, multimapping_hash)
                    read_hash = {}
                    read_counter = 0
                prev_r1name = r1name

            if r1name not in read_hash:
                rec = ReadRecord()
                rec.r1[ensttype] = r1
                rec.flags[ensttype] = enstpriority
                rec.quality = paired_mismatch_score
                rec.mult_ensts[ensttype].append(mapped_enst_full)
                rec.enst[mapped_enst] = mapped_enst_full
                rec.master_enst[ensttype] = mapped_enst_full
                read_hash[r1name] = rec
                continue

            rec = read_hash[r1name]
            if paired_mismatch_score < rec.quality:
                continue
            if paired_mismatch_score > rec.quality:
                newrec = ReadRecord()
                newrec.r1[ensttype] = r1
                newrec.flags[ensttype] = enstpriority
                newrec.quality = paired_mismatch_score
                newrec.mult_ensts[ensttype].append(mapped_enst_full)
                newrec.enst[mapped_enst] = mapped_enst_full
                newrec.master_enst[ensttype] = mapped_enst_full
                read_hash[r1name] = newrec
                continue

            if ensttype in rec.flags:
                if mapped_enst in rec.enst:
                    old_full = rec.enst[mapped_enst]
                    if old_full + '_withgenomeflank' == mapped_enst_full or old_full + '_spliced' == mapped_enst_full:
                        pass
                    elif old_full == mapped_enst_full + '_withgenomeflank' or old_full == mapped_enst_full + '_spliced':
                        rec.r1[ensttype] = r1
                        rec.flags[ensttype] = enstpriority
                        for i, v in enumerate(rec.mult_ensts.get(ensttype, [])):
                            if v == old_full:
                                rec.mult_ensts[ensttype][i] = mapped_enst_full
                        rec.enst[mapped_enst] = mapped_enst_full
                        rec.master_enst[ensttype] = mapped_enst_full
                    else:
                        for i, v in enumerate(rec.mult_ensts.get(ensttype, [])):
                            if v == old_full:
                                rec.mult_ensts[ensttype][i] = mapped_enst_full + '_DOUBLEMAP'
                        rec.master_enst[ensttype] = mapped_enst_full + '_DOUBLEMAP'
                elif enstpriority < rec.flags[ensttype]:
                    rec.r1[ensttype] = r1
                    rec.flags[ensttype] = enstpriority
                    rec.mult_ensts[ensttype].insert(0, mapped_enst_full)
                    rec.enst[mapped_enst] = mapped_enst_full
                    rec.master_enst[ensttype] = mapped_enst_full
                else:
                    rec.mult_ensts[ensttype].append(mapped_enst_full)
            elif ensttype in RRNA_EXTRA_HASH and RRNA_EXTRA_HASH[ensttype] in rec.r1:
                old_rrna = RRNA_EXTRA_HASH[ensttype]
                rec.r1.pop(old_rrna, None)
                rec.flags.pop(old_rrna, None)
                rec.master_enst.pop(old_rrna, None)
                rec.mult_ensts.pop(old_rrna, None)
                rec.enst.pop('NR_046235.1', None)
                rec.r1[ensttype] = r1
                rec.flags[ensttype] = enstpriority
                rec.quality = paired_mismatch_score
                rec.enst[mapped_enst] = mapped_enst_full
                rec.master_enst[ensttype] = mapped_enst_full
                rec.mult_ensts[ensttype].append(mapped_enst_full)
            elif ensttype in RRNA_EXTRA_HASH_REV:
                rna_flag = any(el in rec.r1 for el in RRNA_EXTRA_HASH_REV[ensttype])
                if not rna_flag:
                    rec.r1[ensttype] = r1
                    rec.flags[ensttype] = enstpriority
                    rec.quality = paired_mismatch_score
                    rec.enst[mapped_enst] = mapped_enst_full
                    rec.mult_ensts[ensttype].append(mapped_enst_full)
                    rec.master_enst[ensttype] = mapped_enst_full
            else:
                rec.r1[ensttype] = r1
                rec.flags[ensttype] = enstpriority
                rec.quality = paired_mismatch_score
                rec.enst[mapped_enst] = mapped_enst_full
                rec.mult_ensts[ensttype].append(mapped_enst_full)
                rec.master_enst[ensttype] = mapped_enst_full

        print_output(read_hash, samout, multimapping_hash)

    with multimapping_out.open('w', encoding='utf-8') as mm:
        for key, count in multimapping_hash.items():
            mm.write(f'{key}\t{count}\n')

    done_file.write_text('jobs done\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
