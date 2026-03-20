#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def file_record(path: Path) -> dict:
    h = hashlib.sha1()
    with path.open('rb') as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return {
        'location': path.resolve().as_uri(),
        'basename': path.name,
        'nameroot': path.stem,
        'nameext': ''.join(path.suffixes[-1:]) if path.suffixes else '',
        'class': 'File',
        'checksum': f'sha1${h.hexdigest()}',
        'size': path.stat().st_size,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description='Write legacy-style REPELEMENTMAPPING output JSON')
    ap.add_argument('--out-json', required=True)
    ap.add_argument('--mapping-json', required=True, help='JSON object mapping output key -> file path')
    args = ap.parse_args()

    mapping = json.loads(args.mapping_json)
    out: dict[str, dict] = {}
    for key, value in mapping.items():
        path = Path(value)
        if not path.exists():
            raise FileNotFoundError(f'Missing expected output for {key}: {path}')
        out[key] = file_record(path)

    Path(args.out_json).write_text(json.dumps(out, indent=4) + '\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
