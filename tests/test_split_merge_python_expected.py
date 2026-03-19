from __future__ import annotations

import csv
import gzip
import hashlib
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
PY_SPLIT = REPO / 'bin/python/split_bam_to_subfiles_SEorPE.py'
PY_MERGE = REPO / 'bin/python/merge_multiple_parsed_files.simplified_20191022.py'


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def unzip_to(src_gz: Path, dst: Path) -> None:
    with gzip.open(src_gz, 'rt', encoding='utf-8', errors='replace') as src, dst.open('w', encoding='utf-8') as out:
        out.write(src.read())


def test_split_matches_expected_manifest(tmp_path: Path) -> None:
    sam_src = REPO / 'tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz'
    sam = tmp_path / 'ip.preRmDup.sam.mini.sam'
    unzip_to(sam_src, sam)

    out_dir = tmp_path / 'split'
    out_dir.mkdir()
    subprocess.run(['python3', str(PY_SPLIT), str(sam), 'SE'], cwd=out_dir, check=True)

    manifest = REPO / 'tests/fixtures/mini/expected/split.expected.manifest.tsv'
    with manifest.open('r', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))

    for row in rows:
        p = out_dir / row['filename']
        assert p.exists(), f'missing file: {p.name}'
        assert str(p.stat().st_size) == row['size_bytes']
        assert sha256_file(p) == row['sha256']


def test_merge_matches_expected_output(tmp_path: Path) -> None:
    in1 = REPO / 'tests/fixtures/mini/source/merge_input_1.parsed_v2.txt'
    in2 = REPO / 'tests/fixtures/mini/source/merge_input_2.parsed_v2.txt'
    expected = REPO / 'tests/fixtures/mini/expected/merged.expected.parsed'
    produced = tmp_path / 'merged.parsed'

    subprocess.run(['python3', str(PY_MERGE), str(produced), str(in1), str(in2)], check=True)

    assert produced.stat().st_size == expected.stat().st_size
    assert sha256_file(produced) == sha256_file(expected)
