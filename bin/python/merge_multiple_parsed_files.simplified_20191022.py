#!/usr/bin/env python3
from __future__ import annotations

import argparse
import fcntl
from collections import defaultdict
from pathlib import Path
import re
import sys


READINFO_RE = re.compile(
    r'All\sreads\:\t(\d+)\tPCR\sduplicates\sremoved\:\t(\d+)\tUsable\sRemaining\:\t(\d+)\tUsable\sfrom\sgenomic\smapping\:\t(\d+)\tUsable\sfrom\sfamily\smapping\:\t(\d+)$'
)


def main() -> int:
    ap = argparse.ArgumentParser(description='Merge parsed files')
    ap.add_argument('output_fi')
    ap.add_argument('files', nargs='+')
    args = ap.parse_args()

    output_fi = Path(args.output_fi)
    files = [Path(p) for p in args.files]

    split_fi1 = list(files[0].parts)
    short_fi1 = split_fi1[-1]
    working_dir = Path(*split_fi1[:-1]) if len(split_fi1) > 1 else Path('.')
    failed_jobs_list = working_dir / 'failed_jobs_list.txt'

    print(f'SHORTFI:{short_fi1}', end='')
    print(f'WORKDIR:{working_dir}', end='')

    lock_fi = Path(str(failed_jobs_list) + '.lck')
    lock_fi.parent.mkdir(parents=True, exist_ok=True)
    with lock_fi.open('w') as s:
        fcntl.flock(s.fileno(), fcntl.LOCK_EX)
        with failed_jobs_list.open('a'):
            pass

    read_sums = {
        'all': 0,
        'usable': 0,
        'genomic': 0,
        'repfamily': 0,
        'total': defaultdict(int),
        'element': defaultdict(lambda: {'ensg_all': None, 'ensg_primary': None, 'readnum': 0}),
    }

    for parsed_file in files:
        parse_file(parsed_file, read_sums)

    output_fi.parent.mkdir(parents=True, exist_ok=True)
    with output_fi.open('w', encoding='utf-8') as out:
        usable = read_sums['usable']
        all_reads = read_sums['all']
        genomic = read_sums['genomic']
        repfamily = read_sums['repfamily']

        out.write(f'#READINFO\tAllReads\t{all_reads}\n')
        out.write(f'#READINFO\tUsableReads\t{usable}\t{usable / all_reads if all_reads else 0}\n')
        out.write(f'#READINFO\tGenomicReads\t{genomic}\t{genomic / usable if usable else 0}\n')
        out.write(f'#READINFO\tRepFamilyReads\t{repfamily}\t{repfamily / usable if usable else 0}\n')

        for element, readnum in sorted(read_sums['total'].items(), key=lambda kv: kv[1], reverse=True):
            out.write(f'TOTAL\t{element}\t{readnum}\t{(readnum / usable) if usable else 0}\n')

        items = sorted(
            read_sums['element'].items(),
            key=lambda kv: kv[1]['readnum'],
            reverse=True,
        )
        for element, info in items:
            out.write(
                f"ELEMENT\t{info['ensg_primary']}\t{info['readnum']}\t"
                f"{(info['readnum'] / usable) if usable else 0}\t{element}\t{info['ensg_all']}\n"
            )

    return 0


def parse_file(parsed_fi: Path, read_sums: dict) -> None:
    with parsed_fi.open('r', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            tmp = line.split('\t')
            if not tmp:
                continue
            if tmp[0] == '#READINFO':
                m = READINFO_RE.search(line)
                if m:
                    read_sums['all'] += int(m.group(1))
                elif len(tmp) > 2 and tmp[1] == 'UsableReads':
                    read_sums['usable'] += int(float(tmp[2]))
                elif len(tmp) > 2 and tmp[1] == 'GenomicReads':
                    read_sums['genomic'] += int(float(tmp[2]))
                elif len(tmp) > 2 and tmp[1] == 'RepFamilyReads':
                    read_sums['repfamily'] += int(float(tmp[2]))
                else:
                    print(f"couldn't parse readinfo line {line}", file=sys.stderr)
            elif tmp[0] == '#READINFO2':
                continue
            elif tmp[0] == 'TOTAL':
                _, element, readnum, *_ = tmp
                read_sums['total'][element] += int(float(readnum))
            else:
                if len(tmp) < 6:
                    continue
                _, ensg_primary, readnum, _, enst_all, ensg_all = tmp[:6]
                existing = read_sums['element'][enst_all]
                if existing['ensg_all'] is not None and existing['ensg_all'] != ensg_all:
                    print(
                        f"error - ensg_all mismatch {ensg_all} {existing['ensg_all']}",
                        file=sys.stderr,
                    )
                if existing['ensg_primary'] is not None and existing['ensg_primary'] != ensg_primary:
                    print(
                        f"error - ensg_primary mismatch {ensg_primary} {existing['ensg_primary']}",
                        file=sys.stderr,
                    )
                existing['ensg_all'] = ensg_all
                existing['ensg_primary'] = ensg_primary
                existing['readnum'] += int(float(readnum))


if __name__ == '__main__':
    raise SystemExit(main())
