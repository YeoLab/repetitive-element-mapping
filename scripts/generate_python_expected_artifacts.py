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
MERGE_IN1 = REPO / 'tests/fixtures/mini/source/merge_input_1.parsed_v2.txt'
MERGE_IN2 = REPO / 'tests/fixtures/mini/source/merge_input_2.parsed_v2.txt'
SPLIT_MANIFEST = REPO / 'tests/fixtures/mini/expected/split.expected.manifest.tsv'
MERGED_EXPECTED = REPO / 'tests/fixtures/mini/expected/merged.expected.parsed'


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
    tmp = REPO / 'results/.expected_tmp'
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)

    sam = tmp / 'ip.preRmDup.sam.mini.sam'
    with gzip.open(SOURCE_GZ, 'rt', encoding='utf-8', errors='replace') as src, sam.open('w', encoding='utf-8') as out:
        out.write(src.read())

    split_dir = tmp / 'split'
    split_dir.mkdir()
    subprocess.run(
        ['python3', str(REPO / 'bin/python/split_bam_to_subfiles_SEorPE.py'), str(sam), 'SE'],
        cwd=split_dir,
        check=True,
    )

    SPLIT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with SPLIT_MANIFEST.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['filename', 'size_bytes', 'sha256'])
        for path in sorted(split_dir.glob('*.tmp')):
            writer.writerow([path.name, path.stat().st_size, sha256_file(path)])

    MERGED_EXPECTED.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            'python3',
            str(REPO / 'bin/python/merge_multiple_parsed_files.simplified_20191022.py'),
            str(MERGED_EXPECTED),
            str(MERGE_IN1),
            str(MERGE_IN2),
        ],
        check=True,
    )

    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
