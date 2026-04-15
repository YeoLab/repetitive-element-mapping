from __future__ import annotations

import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
PARSER = REPO / 'bin/python/parse_bowtie2_output_realtime_includemultifamily_PE.py'


def test_parse_pe_runs_on_sam_input(tmp_path: Path) -> None:
    filelist = tmp_path / 'filelist.tsv'
    filelist.write_text(
        'ENST0001|ENST0002\tENSGX\tGENEX\tRNA28S\tx\n'
        'ENST0003\tENSGY\tGENEY\tRNA5S\tx\n',
        encoding='utf-8',
    )

    sam = tmp_path / 'in.sam'
    sam.write_text(
        '@HD\tVN:1.0\n'
        'READ1:UMI\t99\tENST0001\t1\t255\t10M\t=\t1\t0\tAAAAAAAAAA\tIIIIIIIIII\tAS:i:5\n'
        'READ1:UMI\t147\tENST0001\t1\t255\t10M\t=\t1\t0\tTTTTTTTTTT\tIIIIIIIIII\tAS:i:5\n'
        'READ2:UMI\t83\tENST0003\t1\t255\t10M\t=\t1\t0\tCCCCCCCCCC\tIIIIIIIIII\tAS:i:4\n'
        'READ2:UMI\t163\tENST0003\t1\t255\t10M\t=\t1\t0\tGGGGGGGGGG\tIIIIIIIIII\tAS:i:4\n',
        encoding='utf-8',
    )

    out = tmp_path / 'out.sam'
    subprocess.run(
        [
            'python3',
            str(PARSER),
            'dummy_r1.fastq.gz',
            'dummy_r2.fastq.gz',
            'dummy_bowtie_idx',
            str(out),
            str(filelist),
            '--sam-input',
            str(sam),
        ],
        check=True,
    )

    assert out.exists()
    assert (tmp_path / 'out.sam.done').exists()
    assert (tmp_path / 'out.sam.multimapping_deleted').exists()

    body = [l for l in out.read_text(encoding='utf-8').splitlines() if not l.startswith('@')]
    assert body, 'expected emitted mapped lines'
    assert len(body) % 2 == 0
    assert all('ZZ:Z:' in l for l in body)
