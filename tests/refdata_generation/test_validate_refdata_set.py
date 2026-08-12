"""The mouse acceptance suite (issue -475.13 / M-8).

Mouse has no reference outputs, so these checks are all the evidence there is.
That makes the checks themselves load-bearing: each one gets a case proving it
FAILS on the defect it exists to catch, because a check that cannot fail is
indistinguishable from no check at all.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / 'bin' / 'python' / 'refdata_generation' / 'validate_refdata_set.py'

sys.path.insert(0, str(REPO))
from bin.python.refdata_generation.validate_refdata_set import (  # noqa: E402
    REFERENCE_FAMILY_OVERLAP,
    read_bed_sources,
    read_filelist,
)


TRNA_SRC = 'x-tRNAs.fa.list.wgenome_flank.flankNs'


def write_filelist(path, gencode=(), repbase=(), extra_ids=(), mirna=('MI0000060',)):
    """A MASTER_FILELIST in the real block order: Gencode, RepBase, tRNA, miRNA,
    rmsk leftovers."""
    lines = [f'ENST0000{i}.1\tENSG1\tG{i}\t{fam}\tgenelists.{fam}'
             for i, fam in enumerate(gencode)]
    lines += [f'{name}\t{fam}\t{fam}\t{fam}\tSINE' for name, fam in repbase]
    lines += [f'tRNA-Ala-AGC-1-1\ttRNA-Ala-AGC\ttRNA-Ala-AGC\ttRNA\t{TRNA_SRC}']
    for mi in mirna:
        lines += [f'{mi}\tmiRNA\tmiRNA\tmiRNA\tmmu-mir-1',
                  f'{mi}-proximal\tmiRNA\tmiRNA\tmiRNA\tmmu-mir-1-proximal']
    lines += [f'{name}\tX\tX\tX\tLINE' for name in extra_ids]
    path.write_text('\n'.join(lines) + '\n')
    return path


def write_bed(path, rows):
    path.write_text(''.join('\t'.join(map(str, r)) + '\n' for r in rows))
    return path


def write_fasta(path, headers):
    path.write_text(''.join(f'>{h}\nACGT\n' for h in headers))
    return path


@pytest.fixture
def good(tmp_path):
    """A minimal set that passes all six checks."""
    fl = write_filelist(tmp_path / 'fl.tsv',
                        gencode=['RNU1', 'SNORD'],
                        repbase=[('B1_MM', 'SINE'), ('L1MD', 'LINE')],
                        extra_ids=['ALUY'])
    fa = write_fasta(tmp_path / 'idx.fa', ['B1_MM', 'L1MD', 'ENST00000.1'])
    # 2,000 rmsk rows keeps rmsk above the 0.99 floor alongside the other
    # sources, which is the real profile: rmsk is ~99.8% of every BED
    rows = [('chr1', i, i + 50, 'ALUY', 900, '+') for i in range(0, 200000, 100)]
    rows += [('chr1', 10, 20, 'ENST00000.1', '-', '+'),
             ('chr1', 30, 40, 'tRNA-Ala-AGC-1-1', 'tRNA-Ala-AGC', '-'),
             ('chr1', 50, 60, 'MI0000060', 'mmu-mir-1', '+'),
             ('chr1', 40, 50, 'MI0000060-proximal', 'mmu-mir-1-proximal', '+'),
             ('chr1', 60, 70, 'MI0000060-proximal', 'mmu-mir-1-proximal', '+')]
    bed = write_bed(tmp_path / 'u.bed', rows)
    prov = tmp_path / 'prov.tsv'
    prov.write_text('family\tsource_header\tsource_file\tlength\tsha256\n'
                    'B1_MM\tB1_Mm_SINE1/7SL_Mus\tf\t100\tdeadbeef\n'
                    'L1MD\tL1MD_LINE_Mus\tf\t100\tdeadbeef\n')
    return dict(fl=fl, fa=fa, bed=bed, prov=prov, tmp=tmp_path)


def run(good, **override):
    args = [sys.executable, str(SCRIPT), '--assembly', 'test',
            '--master-filelist', str(override.get('fl', good['fl'])),
            '--index-fasta', str(override.get('fa', good['fa'])),
            '--unique-genomic-elements', str(override.get('bed', good['bed']))]
    prov = override.get('prov', good['prov'])
    if prov is not None:
        args += ['--repbase-provenance', str(prov)]
    p = subprocess.run(args, capture_output=True, text=True, cwd=REPO)
    return p.returncode, p.stderr


def test_a_clean_set_passes(good):
    rc, out = run(good)
    assert rc == 0, out
    assert 'PASS (6 passed, 0 failed, 0 skipped)' in out


# --- each check must fail on the defect it exists to catch --------------------

def test_check_a_fails_on_an_index_header_absent_from_the_filelist(good):
    """An index sequence the filelist cannot type aligns reads nothing can label."""
    rc, out = run(good, fa=write_fasta(good['tmp'] / 'bad.fa', ['B1_MM', 'GHOST']))
    assert rc == 1
    assert 'A FAIL' in out


def test_check_b_fails_on_a_bed_name_absent_from_the_filelist(good):
    """T-07: read_peakfi cannot type the peak, so the row is dead weight."""
    rows = [('chr1', 10, 20, 'NOT_IN_FILELIST', 900, '+')]
    rc, out = run(good, bed=write_bed(good['tmp'] / 'bad.bed', rows))
    assert rc == 1
    assert 'B FAIL' in out


def test_check_c_fails_on_a_coordinate_suffix_header(good):
    """T-05: '::' means the repeat portion was cut per genomic instance."""
    rc, out = run(good, fa=write_fasta(good['tmp'] / 'bad.fa',
                                       ['B1_MM::chr1:100-200(+)']))
    assert rc == 1
    assert 'C FAIL' in out


def test_check_d_fails_when_a_family_is_in_both_portions(good):
    """M-5: RNU1 in the RepBase block splits U1 reads with genelists.RNU1."""
    fl = write_filelist(good['tmp'] / 'bad.tsv',
                        gencode=['RNU1'], repbase=[('U1', 'RNU1')], extra_ids=['ALUY'])
    rc, out = run(good, fl=fl, fa=write_fasta(good['tmp'] / 'u.fa', ['U1']),
                  prov=None)
    assert rc == 1
    assert 'D FAIL' in out


def test_check_d_tolerates_exactly_the_overlap_the_hg38_reference_has(good):
    """SNORD and YRNA are in both portions of the reference set itself, so the
    literal 'zero overlap' rule would fail on hg38. Mouse may not ADD to it."""
    assert REFERENCE_FAMILY_OVERLAP == {'SNORD', 'YRNA'}
    fl = write_filelist(good['tmp'] / 'ok.tsv',
                        gencode=['SNORD'], repbase=[('U3', 'SNORD')], extra_ids=['ALUY'])
    _, out = run(good, fl=fl, fa=write_fasta(good['tmp'] / 'u3.fa', ['U3']), prov=None)
    assert 'D PASS' in out


def test_check_e_fails_on_simple_repeat_rows(good):
    """T-07 rule 1: the trf track is not a source; the reference has 0 rows."""
    rows = [('chr1', i, i + 50, 'ALUY', 900, '+') for i in range(0, 200000, 100)]
    rows += [('chr1', 1, 2, 'trf', 102, '.')]
    rc, out = run(good, bed=write_bed(good['tmp'] / 'bad.bed', rows))
    assert rc == 1
    assert 'E FAIL' in out


def test_check_e_fails_when_an_optional_source_is_silently_missing(good):
    """Every optional source is a flag the generator only WARNS about. mm39's
    superseded BED had zero tRNA and zero miRNA rows because neither input
    existed, and nothing failed."""
    rows = [('chr1', i, i + 50, 'ALUY', 900, '+') for i in range(0, 200000, 100)]
    rows += [('chr1', 10, 20, 'ENST00000.1', '-', '+')]
    rc, out = run(good, bed=write_bed(good['tmp'] / 'bad.bed', rows))
    assert rc == 1
    assert 'E FAIL' in out
    assert 'no trna rows at all' in out and 'no mirna rows at all' in out


def test_check_f_fails_when_a_repeat_family_has_no_provenance_row(good):
    """The other side of the T-05 guard: an emitted family must trace to the
    RepBase record it came from."""
    prov = good['tmp'] / 'short.tsv'
    prov.write_text('family\tsource_header\tsource_file\tlength\tsha256\n'
                    'B1_MM\tB1_Mm_SINE1/7SL_Mus\tf\t100\tdeadbeef\n')
    rc, out = run(good, prov=prov)
    assert rc == 1
    assert 'F FAIL' in out


def test_check_f_is_skipped_not_failed_without_a_sidecar(good):
    """The hg38 reference index is a downloaded 2020 artifact with no sidecar,
    and the suite has to pass on it."""
    rc, out = run(good, prov=None)
    assert rc == 0
    assert 'F SKIP' in out


# --- block segmentation, which every family-level check depends on ------------

def test_repbase_block_is_bounded_by_the_trna_source_list_not_column_4(tmp_path):
    """mm10's MamSINE1 row carries family 'tRNA' (mm39's carries 'tRNA-RTE'),
    straight from each assembly's rmsk table. Delimiting on column 4 cut mm10's
    RepBase block at 387 rows instead of 1,071."""
    fl = write_filelist(tmp_path / 'fl.tsv', gencode=['RNU1'],
                        repbase=[('MAMSINE1', 'tRNA'), ('B1_MM', 'SINE')])
    _, _, repbase_families, repbase_names = read_filelist(fl)
    assert repbase_names == {'MAMSINE1', 'B1_MM'}
    assert repbase_families == {'tRNA', 'SINE'}


def test_pipe_joined_filelist_ids_are_split(tmp_path):
    """Column 1 holds pipe-joined ids; read_in_filelists splits them, so an
    index header matching one component must count as resolved."""
    fl = tmp_path / 'fl.tsv'
    fl.write_text('ENST1.1|ENST2.1\tENSG1\tG\tRNU1\tgenelists.RNU1\n'
                  f'B1_MM\tSINE\tSINE\tSINE\tSINE\n'
                  f'tRNA-Ala-AGC-1-1\ta\ta\ttRNA\t{TRNA_SRC}\n')
    ids, _, _, _ = read_filelist(fl)
    assert {'ENST1.1', 'ENST2.1'} <= ids


def test_bed_source_classification_matches_the_generator(tmp_path):
    bed = write_bed(tmp_path / 'b.bed', [
        ('chr1', 1, 2, 'ALUY', 900, '+'),                        # rmsk
        ('chr1', 1, 2, 'ENST1.1', '-', '+'),                     # gencode
        ('chr1', 1, 2, 'tRNA-Ala-AGC-1-1', 'tRNA-Ala-AGC', '-'),  # trna
        ('chr1', 1, 2, 'MI0000060', 'mmu-mir-1', '+'),           # mirna
        ('chr1', 1, 2, 'trf', 102, '.'),                         # simple repeat
    ])
    total, unresolved, counts, mirna_ids = read_bed_sources(bed, set())
    assert total == 5 and unresolved == 5
    assert dict(counts) == {'rmsk': 1, 'gencode': 1, 'trna': 1,
                            'mirna': 1, 'simplerepeat': 1}
    assert mirna_ids == {'MI0000060'}
