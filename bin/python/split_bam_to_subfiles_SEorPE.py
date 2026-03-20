#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


PREFIXES = [a + b for a in 'ACGTN' for b in 'ACGTN']


def sam_stream(path: Path):
    if path.suffix == '.sam':
        with path.open('r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                yield line
    elif path.suffix == '.bam':
        samtools_bin = os.environ.get('SAMTOOLS_BIN') or shutil.which('samtools')
        if not samtools_bin:
            repo_samtools = Path(__file__).resolve().parents[2] / '.conda-env' / 'bin' / 'samtools'
            samtools_bin = str(repo_samtools) if repo_samtools.exists() else 'samtools'
        proc = subprocess.Popen(
            [samtools_bin, 'view', '-h', str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line
        stderr = proc.stderr.read() if proc.stderr else ''
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError(f'samtools view failed ({ret}): {stderr}')
    else:
        raise ValueError(f'weird - {path} not either sam or bam file format - exit')


def main() -> int:
    ap = argparse.ArgumentParser(description='Split SAM/BAM into 25 UMI-prefix tmp files')
    ap.add_argument('sam_fi')
    ap.add_argument('se_or_pe_flag', choices=['SE', 'PE'])
    args = ap.parse_args()

    sam_path = Path(args.sam_fi)
    sam_short = sam_path.name
    mode = args.se_or_pe_flag

    if mode == 'PE':
        print('splitting in paired-end mode', file=sys.stderr)
    else:
        print('splitting in single-end mode', file=sys.stderr)

    filehandles = {
        p: (Path(f'{p}.{sam_short}.tmp')).open('w', encoding='utf-8')
        for p in PREFIXES
    }

    try:
        stream = sam_stream(sam_path)
        it = iter(stream)
        for r1 in it:
            if r1.startswith('@'):
                continue
            r1 = r1.rstrip('\n')
            tmp_r1 = r1.split('\t')
            r1_name = tmp_r1[0].split()[0]
            r1_flag = int(tmp_r1[1])

            r2 = None
            tmp_r2: list[str] | None = None
            if mode == 'PE':
                try:
                    r2 = next(it).rstrip('\n')
                except StopIteration:
                    break
                tmp_r2 = r2.split('\t')
                r2_name = tmp_r2[0].split()[0]
                if r1_name != r2_name:
                    print(
                        f'paired end mismatch error: {sam_path} r1 {tmp_r1[0]} r2 {tmp_r2[0]}',
                        file=sys.stderr,
                    )
                if not r1_flag:
                    print(f'error {r1} {r2}', file=sys.stderr)

            if mode == 'PE':
                if r1_flag in (77, 141):
                    continue
            else:
                if r1_flag == 4:
                    continue

            if mode == 'PE':
                if r1_flag in (99, 355, 83, 339):
                    pass
                elif r1_flag in (147, 403, 163, 419):
                    assert r2 is not None and tmp_r2 is not None
                    tmp_r1, tmp_r2 = r2.split('\t'), r1.split('\t')
                else:
                    continue

                randommer = tmp_r1[0].split(':')[0]
                first2 = randommer[:2]
                if first2 in filehandles and r2 is not None:
                    filehandles[first2].write(f'{r1}\n{r2}\n')
            else:
                if r1_flag not in (16, 272, 0, 256):
                    continue
                read_name_parts = tmp_r1[0].split('_')
                randommer = read_name_parts[-1]
                first2 = randommer[:2]
                if first2 in filehandles:
                    filehandles[first2].write(f'{r1}\n')
    finally:
        for fh in filehandles.values():
            fh.close()

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
