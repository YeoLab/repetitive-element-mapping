import math
import subprocess
from pathlib import Path


def parse_merged(path: Path):
    headers = {}
    totals = {}
    elements = {}

    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            cols = raw.rstrip('\n').split('\t')
            if cols[0] == '#READINFO':
                key = cols[1]
                count = int(float(cols[2]))
                frac = float(cols[3]) if len(cols) > 3 else None
                headers[key] = (count, frac)
            elif cols[0] == 'TOTAL':
                totals[cols[1]] = int(float(cols[2]))
            elif cols[0] == 'ELEMENT':
                elements[cols[4]] = (cols[1], int(float(cols[2])), cols[5])
    return headers, totals, elements


py_out = Path(snakemake.input['py'])
in1 = Path(snakemake.input['in1'])
in2 = Path(snakemake.input['in2'])
out_ok = Path(snakemake.output[0])

perl_out = py_out.parent / 'merged.perl.parsed'
perl_script = Path(snakemake.params['perl_script'])
subprocess.run(['perl', str(perl_script), str(perl_out), str(in1), str(in2)], check=True)

perl_h, perl_t, perl_e = parse_merged(perl_out)
py_h, py_t, py_e = parse_merged(py_out)

assert perl_t == py_t
assert perl_e == py_e
assert perl_h.keys() == py_h.keys()
for key in perl_h:
    p_count, p_frac = perl_h[key]
    y_count, y_frac = py_h[key]
    assert p_count == y_count
    if p_frac is not None and y_frac is not None:
        assert math.isclose(p_frac, y_frac, rel_tol=1e-12, abs_tol=1e-12)

out_ok.parent.mkdir(parents=True, exist_ok=True)
out_ok.write_text('ok\n', encoding='utf-8')
