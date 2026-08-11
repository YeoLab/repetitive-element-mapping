"""Selection-rule tests for generate_unique_genomic_elements.py (issue -475.2).

Every source in the UniqueGenomicElements BED gets an inclusion and an exclusion case.
Rules were derived by per-source diff against the hg38 reference BED; see
docs/CHANGELOG-refdata-validation-2026-07-25.md §T-07.
"""

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / 'bin' / 'python' / 'refdata_generation' / 'generate_unique_genomic_elements.py'

sys.path.insert(0, str(REPO))
from bin.python.refdata_generation.generate_unique_genomic_elements import (  # noqa: E402
    GENCODE_EXCLUDED_FAMILIES,
    parse_trna_bed_rows,
    read_chrom_allowlist,
    read_gencode_allowed_transcripts,
)


def gtf_line(chrom, start, end, score, strand, gene_id, transcript_id=None):
    attrs = f'gene_id "{gene_id}"; transcript_id "{transcript_id or gene_id}";'
    return f'{chrom}\tsrc\texon\t{start}\t{end}\t{score}\t{strand}\t.\t{attrs}\n'


@pytest.fixture
def inputs(tmp_path):
    """Minimal one-row-per-case inputs exercising both sides of every rule."""
    rm = tmp_path / 'rm.tsv.gz'
    with gzip.open(rm, 'wt') as fh:
        fh.write(gtf_line('chr1', 101, 200, '1892.000000', '+', 'L1P5'))        # kept
        fh.write(gtf_line('chr1_ML143354v1_fix', 1, 50, '900.0', '-', 'AluY'))  # dropped: scaffold

    trna = tmp_path / 'trna.tsv.gz'
    with gzip.open(trna, 'wt') as fh:
        fh.write(gtf_line('chr1', 301, 380, '1000.000000', '-', 'tRNA-Asn-GTT-2-3'))
        fh.write(gtf_line('chrUn_novel', 1, 70, '1000.0', '+', 'tRNA-Val-TAC-3-1'))  # dropped

    ucsc = tmp_path / 'parsed_ucsc.tsv'
    cols = lambda tid, start, end: '\t'.join(  # noqa: E731
        ['ENSG00000000001.1', tid, 'chr1', '+', str(start), str(end), '0', '0', '1',
         f'{start},', f'{end},']) + '\n'
    ucsc.write_text(
        cols('ENST00000383925.1', 500, 600) +   # in MASTER_FILELIST, family RNU1 -> kept
        cols('ENST00000600000.1', 700, 800) +   # in MASTER_FILELIST, family RNA5S -> excluded
        cols('ENST00000999999.1', 900, 999)     # absent from MASTER_FILELIST -> excluded
    )

    master = tmp_path / 'master.tsv'
    master.write_text(
        'ENST00000383925.1\tENSG00000206652.1\tRNU1-1\tRNU1\tgenelists.RNU1\n'
        'ENST00000600000.1\tENSG00000000002.1\tRNA5S1\tRNA5S\tgenelists.RNA5S\n'
    )

    gff3 = tmp_path / 'mirna.gff3'
    gff3.write_text(
        '# comment\n'
        'chr1\t.\tmiRNA_primary_transcript\t2001\t2060\t.\t+\t.\tID=MI0000060;Name=hsa-mir-1\n'
        'chrUn_novel\t.\tmiRNA_primary_transcript\t10\t70\t.\t+\t.\tID=MI0009999;Name=hsa-mir-x\n'
    )

    allowlist = tmp_path / 'chroms.txt'
    allowlist.write_text('chr1\n')

    return dict(rm=rm, trna=trna, ucsc=ucsc, master=master, gff3=gff3,
                allowlist=allowlist, out=tmp_path / 'out.bed')


def run(inputs, tmp_path):
    subprocess.run([sys.executable, str(SCRIPT),
                    '--repeatmasker', str(inputs['rm']),
                    '--trna', str(inputs['trna']),
                    '--gff3', str(inputs['gff3']),
                    '--parsed-ucsc', str(inputs['ucsc']),
                    '--master-filelist', str(inputs['master']),
                    '--chrom-allowlist', str(inputs['allowlist']),
                    '--assembly', 'test',
                    '--output', str(inputs['out'])],
                   check=True, cwd=REPO, capture_output=True)
    return [ln.split('\t') for ln in inputs['out'].read_text().splitlines()]


# --- simple repeats: the source must not be accepted at all --------------------------

def test_simplerepeats_flag_is_rejected(inputs, tmp_path):
    """The hg38 reference contains zero simple-repeat rows; the option must be gone."""
    proc = subprocess.run([sys.executable, str(SCRIPT),
                           '--repeatmasker', str(inputs['rm']),
                           '--master-filelist', str(inputs['master']),
                           '--chrom-allowlist', str(inputs['allowlist']),
                           '--assembly', 'test',
                           '--output', str(inputs['out']),
                           '--simplerepeats', 'x'],
                          capture_output=True, text=True, cwd=REPO)
    assert proc.returncode != 0
    assert 'unrecognized arguments: --simplerepeats' in proc.stderr


# --- RepeatMasker -------------------------------------------------------------------

def test_repeatmasker_row_included_with_numeric_score(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    hits = [r for r in rows if r[3] == 'L1P5']
    assert hits == [['chr1', '100', '200', 'L1P5', '1892', '+']]


def test_repeatmasker_row_on_unlisted_scaffold_excluded(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    assert not [r for r in rows if r[0] == 'chr1_ML143354v1_fix']


# --- tRNA ---------------------------------------------------------------------------

def test_trna_col5_is_family_not_score(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    hits = [r for r in rows if r[3] == 'tRNA-Asn-GTT-2-3']
    assert hits == [['chr1', '300', '380', 'tRNA-Asn-GTT-2-3', 'tRNA-Asn-GTT', '-']]


def test_trna_row_on_unlisted_scaffold_excluded(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    assert not [r for r in rows if r[3] == 'tRNA-Val-TAC-3-1']


@pytest.mark.parametrize('gene_id,family', [
    ('tRNA-Asn-GTT-2-3', 'tRNA-Asn-GTT'),
    ('tRNA-Val-TAC-3-1', 'tRNA-Val-TAC'),
    ('nmt-tRNA-Gln-CTG-1-1', 'nmt-tRNA-Gln-CTG'),
    ('tRNA-Met', 'tRNA-Met'),                       # no copy suffix -> unchanged
    ('nm-tRNA-Tyr-GTA-chr10-3', 'nm-tRNA-Tyr-GTA-chr10-3'),  # non-numeric tail -> unchanged
])
def test_trna_family_derivation(tmp_path, gene_id, family):
    src = tmp_path / 't.tsv'
    src.write_text(gtf_line('chr1', 1, 10, '1000.0', '+', gene_id))
    assert list(parse_trna_bed_rows(src))[0][4] == family


# --- Gencode ------------------------------------------------------------------------

def test_gencode_transcript_in_master_filelist_included(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    hits = [r for r in rows if r[3] == 'ENST00000383925.1']
    assert hits == [['chr1', '500', '600', 'ENST00000383925.1', '-', '+']]


def test_gencode_transcript_absent_from_master_filelist_excluded(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    assert not [r for r in rows if r[3] == 'ENST00000999999.1']


def test_gencode_rrna_family_excluded(inputs, tmp_path):
    """RNA5S/RNA5-8S/MT* are multi-copy; they belong to repeat families, not the
    unique genome. 582 of the 591 hg38 omissions are these."""
    rows = run(inputs, tmp_path)
    assert not [r for r in rows if r[3] == 'ENST00000600000.1']


def test_gencode_allowlist_drops_every_excluded_family(tmp_path):
    master = tmp_path / 'm.tsv'
    master.write_text(''.join(
        f'ENST0000000{i}.1\tENSG1\tG{i}\t{fam}\tgenelists.{fam}\n'
        for i, fam in enumerate(sorted(GENCODE_EXCLUDED_FAMILIES) + ['RNU1'])
    ))
    allowed = read_gencode_allowed_transcripts(master)
    assert allowed == {f'ENST0000000{len(GENCODE_EXCLUDED_FAMILIES)}.1'}


# --- miRNA --------------------------------------------------------------------------

def test_mirna_emits_primary_plus_two_proximal(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    hits = sorted(r for r in rows if r[3].startswith('MI0000060'))
    assert hits == [
        ['chr1', '1500', '2000', 'MI0000060-proximal', 'hsa-mir-1-proximal', '+'],
        ['chr1', '2000', '2060', 'MI0000060', 'hsa-mir-1', '+'],
        ['chr1', '2060', '2560', 'MI0000060-proximal', 'hsa-mir-1-proximal', '+'],
    ]


def test_mirna_on_unlisted_scaffold_excluded(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    assert not [r for r in rows if r[3].startswith('MI0009999')]


# --- allowlist parsing --------------------------------------------------------------

def test_chrom_allowlist_accepts_fai(tmp_path):
    fai = tmp_path / 'g.fa.fai'
    fai.write_text('chr1\t248956422\t112\t60\t61\nchr2\t242193529\t1\t60\t61\n')
    assert read_chrom_allowlist(fai) == {'chr1', 'chr2'}


def test_chrom_allowlist_skips_blanks_and_comments(tmp_path):
    lst = tmp_path / 'c.txt'
    lst.write_text('# scaffolds\nchr1\n\nchr2\n')
    assert read_chrom_allowlist(lst) == {'chr1', 'chr2'}


# --- schema -------------------------------------------------------------------------

def test_output_is_six_column_bed_with_nonnegative_starts(inputs, tmp_path):
    rows = run(inputs, tmp_path)
    assert rows
    assert all(len(r) == 6 for r in rows)
    assert all(int(r[1]) >= 0 and int(r[2]) > int(r[1]) for r in rows)


# ── miRBase seqid normalization (-475.6) ─────────────────────────────────

@pytest.mark.parametrize("seqid, expected", [
    ("1", "chr1"),
    ("X", "chrX"),
    ("MT", "chrM"),
    ("chr2", "chr2"),        # already UCSC-named; left alone
    ("chrX", "chrX"),
    ("GL456210.1", "chrGL456210.1"),
])
def test_ucsc_chrom_normalizes_gff3_seqids(seqid, expected):
    from generate_unique_genomic_elements import ucsc_chrom
    assert ucsc_chrom(seqid) == expected


def test_gff3_mirna_seqids_are_normalized(tmp_path):
    """miRBase v23 mmu.gff3 mixes conventions: 1,164 bare rows and 26 'chr' rows.

    Un-normalized seqids silently produce intervals matching no chromosome,
    because every other artifact here is UCSC-named.
    """
    from generate_unique_genomic_elements import parse_gff3_mirna
    f = tmp_path / "mmu.gff3"
    f.write_text(
        "##gff-version 3\n"
        "1\t.\tmiRNA_primary_transcript\t10\t20\t.\t+\t.\tID=MI1;Name=mmu-mir-1\n"
        "chr2\t.\tmiRNA_primary_transcript\t30\t40\t.\t-\t.\tID=MI2;Name=mmu-mir-2\n"
    )
    assert [r[0] for r in parse_gff3_mirna(f)] == ["chr1", "chr2"]
