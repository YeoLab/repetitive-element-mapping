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


produced = Path(snakemake.input['produced'])
expected = Path(snakemake.input['expected'])
out_ok = Path(snakemake.output[0])

assert produced.exists(), f'missing produced file: {produced}'
assert expected.exists(), f'missing expected file: {expected}'
assert produced.stat().st_size == expected.stat().st_size, 'size mismatch'
assert sha256_file(produced) == sha256_file(expected), 'sha mismatch'

out_ok.parent.mkdir(parents=True, exist_ok=True)
out_ok.write_text('ok\n', encoding='utf-8')
