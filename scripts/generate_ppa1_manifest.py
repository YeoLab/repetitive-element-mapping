#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        while True:
            data = fh.read(chunk_size)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description='Generate expected output manifest for PPA1 example files')
    parser.add_argument('--input-dir', default='data/PPA1_rep1/results', help='Directory containing expected result files')
    parser.add_argument('--output', default='docs/ppa1_expected_outputs_manifest.tsv', help='TSV output path')
    args = parser.parse_args()

    in_dir = Path(args.input_dir)
    out = Path(args.output)

    files = sorted(p for p in in_dir.iterdir() if p.is_file())
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open('w', encoding='utf-8') as fh:
        fh.write('relative_path\tsize_bytes\tsha256\n')
        for path in files:
            fh.write(f"{path.as_posix()}\t{path.stat().st_size}\t{sha256_file(path)}\n")


if __name__ == '__main__':
    main()
