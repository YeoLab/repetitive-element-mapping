from __future__ import annotations

import csv
import hashlib
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def test_mini_fixture_manifest_matches_files() -> None:
    manifest = Path('tests/fixtures/mini/manifest.tsv')
    assert manifest.exists(), 'manifest missing; run scripts/make_mini_fixtures.py'

    with manifest.open('r', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = list(reader)

    assert rows, 'manifest contains no fixture rows'

    for row in rows:
        path = Path(row['relative_path'])
        assert path.exists(), f'missing fixture file: {path}'
        assert str(path.stat().st_size) == row['size_bytes'], f'size mismatch: {path}'
        assert sha256_file(path) == row['sha256'], f'sha256 mismatch: {path}'
