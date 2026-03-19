#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import shutil
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SOURCE_GZ = REPO / 'tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz'
MANIFEST = REPO / 'tests/fixtures/mini/expected/split.perl.manifest.tsv'


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    tmp = REPO / 'results/.split_manifest_tmp'
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)

    sam = tmp / 'ip.preRmDup.sam.mini.sam'
    with gzip.open(SOURCE_GZ, 'rt', encoding='utf-8', errors='replace') as src, sam.open('w', encoding='utf-8') as out:
        out.write(src.read())

    subprocess.run(['perl', str(REPO / 'bin/perl/split_bam_to_subfiles_SEorPE.pl'), str(sam), 'SE'], cwd=tmp, check=True)

    files = sorted(tmp.glob('*.tmp'))
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['filename', 'size_bytes', 'sha256'])
        for path in files:
            writer.writerow([path.name, path.stat().st_size, sha256_file(path)])

    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
