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


split_dir = Path(snakemake.input['split_dir'])
manifest = Path(snakemake.input['manifest'])
out_ok = Path(snakemake.output[0])

with manifest.open('r', encoding='utf-8') as fh:
    rows = list(csv.DictReader(fh, delimiter='\t'))

for row in rows:
    p = split_dir / row['filename']
    assert p.exists(), f'missing split output {p}'
    assert str(p.stat().st_size) == row['size_bytes'], f'size mismatch {p}'
    assert sha256_file(p) == row['sha256'], f'sha mismatch {p}'

out_ok.parent.mkdir(parents=True, exist_ok=True)
out_ok.write_text('ok\n', encoding='utf-8')
