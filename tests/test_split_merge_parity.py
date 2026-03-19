from __future__ import annotations

import gzip
import hashlib
import math
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
PERL_SPLIT = REPO / 'bin/perl/split_bam_to_subfiles_SEorPE.pl'
PY_SPLIT = REPO / 'bin/python/split_bam_to_subfiles_SEorPE.py'
PERL_MERGE = REPO / 'bin/perl/merge_multiple_parsed_files.simplified_20191022.pl'
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


def test_split_script_parity(tmp_path: Path) -> None:
    sam_src = REPO / 'tests/fixtures/mini/source/ip.preRmDup.sam.mini.gz'
    sam = tmp_path / 'ip.preRmDup.sam'
    unzip_to(sam_src, sam)

    perl_dir = tmp_path / 'perl_split'
    py_dir = tmp_path / 'py_split'
    perl_dir.mkdir()
    py_dir.mkdir()

    subprocess.run(['perl', str(PERL_SPLIT), str(sam), 'SE'], cwd=perl_dir, check=True)
    subprocess.run(['python3', str(PY_SPLIT), str(sam), 'SE'], cwd=py_dir, check=True)

    perl_files = sorted(perl_dir.glob('*.tmp'))
    py_files = sorted(py_dir.glob('*.tmp'))

    assert len(perl_files) == 25
    assert len(py_files) == 25
    assert [p.name for p in perl_files] == [p.name for p in py_files]

    for p_perl, p_py in zip(perl_files, py_files):
        assert sha256_file(p_perl) == sha256_file(p_py), f'mismatch: {p_perl.name}'


def parse_merged(path: Path):
    headers: dict[str, tuple[int, float | None]] = {}
    totals: dict[str, int] = {}
    elements: dict[str, tuple[str, int, str]] = {}

    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            cols = line.split('\t')
            if cols[0] == '#READINFO':
                key = cols[1]
                count = int(float(cols[2]))
                frac = float(cols[3]) if len(cols) > 3 else None
                headers[key] = (count, frac)
            elif cols[0] == 'TOTAL':
                totals[cols[1]] = int(float(cols[2]))
            elif cols[0] == 'ELEMENT':
                # ELEMENT primary readnum frac element_id gene_id
                elements[cols[4]] = (cols[1], int(float(cols[2])), cols[5])

    return headers, totals, elements


def test_merge_script_parity(tmp_path: Path) -> None:
    ip_parsed = tmp_path / 'ip.prefix.parsed_v2.20201210.txt'
    input_parsed = tmp_path / 'input.prefix.parsed_v2.20201210.txt'

    ip_parsed.write_text(
        '\n'.join(
            [
                '#READINFO\tAll reads:\t100\tPCR duplicates removed:\t20\tUsable Remaining:\t80\tUsable from genomic mapping:\t30\tUsable from family mapping:\t50',
                '#READINFO\tUsableReads\t80',
                '#READINFO\tGenomicReads\t30\t0.375',
                '#READINFO\tRepFamilyReads\t50\t0.625',
                'TOTAL\tRNA18S\t20\t250000',
                'TOTAL\tRNA28S\t30\t375000',
                'ELEMENT\tRNA18S\t12\t150000\tRNA18S||ENST1\tGENE1',
                'ELEMENT\tRNA28S\t18\t225000\tRNA28S||ENST2\tGENE2',
            ]
        )
        + '\n',
        encoding='utf-8',
    )
    input_parsed.write_text(
        '\n'.join(
            [
                '#READINFO\tAll reads:\t200\tPCR duplicates removed:\t40\tUsable Remaining:\t160\tUsable from genomic mapping:\t60\tUsable from family mapping:\t100',
                '#READINFO\tUsableReads\t160',
                '#READINFO\tGenomicReads\t60\t0.375',
                '#READINFO\tRepFamilyReads\t100\t0.625',
                'TOTAL\tRNA18S\t40\t250000',
                'TOTAL\tRNA28S\t60\t375000',
                'ELEMENT\tRNA18S\t25\t156250\tRNA18S||ENST1\tGENE1',
                'ELEMENT\tRNA28S\t35\t218750\tRNA28S||ENST2\tGENE2',
            ]
        )
        + '\n',
        encoding='utf-8',
    )

    perl_out = tmp_path / 'merged.perl.parsed'
    py_out = tmp_path / 'merged.py.parsed'

    subprocess.run(
        ['perl', str(PERL_MERGE), str(perl_out), str(ip_parsed), str(input_parsed)],
        check=True,
    )
    subprocess.run(
        ['python3', str(PY_MERGE), str(py_out), str(ip_parsed), str(input_parsed)],
        check=True,
    )

    perl_h, perl_t, perl_e = parse_merged(perl_out)
    py_h, py_t, py_e = parse_merged(py_out)

    assert perl_t == py_t
    assert perl_e == py_e

    assert perl_h.keys() == py_h.keys()
    for key in perl_h:
        p_count, p_frac = perl_h[key]
        y_count, y_frac = py_h[key]
        assert p_count == y_count
        if p_frac is None or y_frac is None:
            assert p_frac == y_frac
        else:
            assert math.isclose(p_frac, y_frac, rel_tol=1e-12, abs_tol=1e-12), key
