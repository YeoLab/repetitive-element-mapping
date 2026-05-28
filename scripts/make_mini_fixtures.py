#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path
from typing import Iterable


def iter_lines_gz(path: Path) -> Iterable[str]:
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            yield line


def write_gz(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8') as out:
        for line in lines:
            out.write(line)


def head_lines(lines: Iterable[str], n: int) -> list[str]:
    out: list[str] = []
    for line in lines:
        out.append(line)
        if len(out) >= n:
            break
    return out


def parsed_subset(lines: Iterable[str], total_rows: int, element_rows: int) -> list[str]:
    out: list[str] = []
    total_seen = 0
    element_seen = 0
    for line in lines:
        if line.startswith('#READINFO'):
            out.append(line)
            continue
        if line.startswith('TOTAL\t') and total_seen < total_rows:
            out.append(line)
            total_seen += 1
            continue
        if line.startswith('ELEMENT\t') and element_seen < element_rows:
            out.append(line)
            element_seen += 1
            continue
        if total_seen >= total_rows and element_seen >= element_rows:
            break
    return out


def tsv_subset(lines: Iterable[str], n_data_rows: int) -> list[str]:
    out: list[str] = []
    for line in lines:
        out.append(line)
        if len(out) == 1:
            continue
        if len(out) - 1 >= n_data_rows:
            break
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        while True:
            data = fh.read(1024 * 1024)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser(description='Create mini deterministic fixtures from PPA1 example outputs')
    p.add_argument('--source-dir', default='data/PPA1_rep1/results')
    p.add_argument('--out-dir', default='tests/fixtures/mini')
    p.add_argument('--sam-lines', type=int, default=3000)
    p.add_argument('--parsed-total-rows', type=int, default=60)
    p.add_argument('--parsed-element-rows', type=int, default=200)
    p.add_argument('--tsv-rows', type=int, default=300)
    args = p.parse_args()

    source = Path(args.source_dir)
    out = Path(args.out_dir)
    src_out = out / 'source'
    expected_out = out / 'expected'

    mapping = {
        'MP2.PPA1_IP1.umi.r1.fqTrTr.sorted.fq.barcode1.preRmDup.sam.gz': ('source/ip.preRmDup.sam.mini.gz', 'sam'),
        'MP2.PPA1_IN1.umi.r1.fqTrTr.sorted.fq.input.preRmDup.sam.gz': ('source/input.preRmDup.sam.mini.gz', 'sam'),
        'MP2.PPA1_IP1.umi.r1.fqTrTr.sorted.fq.barcode1.rmDup.sam.gz': ('expected/ip.rmDup.sam.mini.gz', 'sam'),
        'MP2.PPA1_IN1.umi.r1.fqTrTr.sorted.fq.input.rmDup.sam.gz': ('expected/input.rmDup.sam.mini.gz', 'sam'),
        'MP2.PPA1_IP1.umi.r1.fqTrTr.sorted.fq.barcode1.parsed.gz': ('expected/ip.parsed.mini.gz', 'parsed'),
        'MP2.PPA1_IN1.umi.r1.fqTrTr.sorted.fq.input.parsed.gz': ('expected/input.parsed.mini.gz', 'parsed'),
        'MP2.PPA1_IP1.umi.r1.fqTrTr.sorted.fq.barcode1.rmDup.sam.reparsed.nopipes.tsv.gz': ('expected/ip.reparsed.nopipes.mini.tsv.gz', 'tsv'),
        'MP2.PPA1_IP1.umi.r1.fqTrTr.sorted.fq.barcode1.rmDup.sam.reparsed.withpipes.tsv.gz': ('expected/ip.reparsed.withpipes.mini.tsv.gz', 'tsv'),
    }

    generated: list[Path] = []
    for src_name, (rel_dst, kind) in mapping.items():
        src = source / src_name
        dst = out / rel_dst
        lines = iter_lines_gz(src)
        if kind == 'sam':
            subset = head_lines(lines, args.sam_lines)
        elif kind == 'parsed':
            subset = parsed_subset(lines, args.parsed_total_rows, args.parsed_element_rows)
        elif kind == 'tsv':
            subset = tsv_subset(lines, args.tsv_rows)
        else:
            raise ValueError(kind)
        write_gz(dst, subset)
        generated.append(dst)

    manifest = out / 'manifest.tsv'
    with manifest.open('w', encoding='utf-8') as fh:
        fh.write('relative_path\tsize_bytes\tsha256\n')
        for path in sorted(generated):
            fh.write(f"{path.as_posix()}\t{path.stat().st_size}\t{sha256_file(path)}\n")


if __name__ == '__main__':
    main()
