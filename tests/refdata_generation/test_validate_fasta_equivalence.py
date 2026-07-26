"""Tests for validate_fasta_equivalence.py (issue -475.5)."""

import hashlib
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from bin.python.refdata_generation.validate_fasta_equivalence import (  # noqa: E402
    canonical_digest,
    compare,
    load,
    main,
)


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# --- canonicalization: uppercase BEFORE ambiguity->N ---------------------------------

def test_lowercase_iupac_becomes_n_not_uppercase_iupac():
    """The regression this validator exists for: the old AWK sketch substituted
    [RYMKSWBDHV] before uppercasing, so 'rymk' survived as 'RYMK'."""
    assert canonical_digest([b'rymk']) == canonical_digest([b'NNNN'])
    assert canonical_digest([b'rymk']) == canonical_digest([b'RYMK'])
    # and the buggy behaviour is genuinely distinguishable: RYMK-left-as-is != NNNN
    assert canonical_digest([b'rymk']) != hashlib.sha256(b'RYMK').hexdigest()


@pytest.mark.parametrize('seq', [b'RYMKSWBDHV', b'rymkswbdhv', b'RyMkSwBdHv'])
def test_all_iupac_codes_map_to_n(seq):
    assert canonical_digest([seq]) == canonical_digest([b'N' * 10])


def test_acgt_preserved_and_case_folded():
    assert canonical_digest([b'acgt']) == canonical_digest([b'ACGT'])
    assert canonical_digest([b'ACGT']) != canonical_digest([b'ACGA'])


def test_multiline_records_hash_same_as_single_line(tmp_path):
    multi = write(tmp_path, 'a.fa', '>x\nACGT\nacgt\nRY\n')
    single = write(tmp_path, 'b.fa', '>x\nACGTACGTNN\n')
    assert load(multi)['X'] == load(single)['X']


def test_lowercase_multiline_fasta_matches_uppercase(tmp_path):
    lower = write(tmp_path, 'l.fa', '>Mer5a\nacgtry\nmkacgt\n')
    upper = write(tmp_path, 'u.fa', '>MER5A\nACGTNNNNACGT\n')
    assert load(lower) == load(upper)


# --- duplicate headers ---------------------------------------------------------------

def test_duplicate_headers_raise(tmp_path):
    dup = write(tmp_path, 'd.fa', '>x\nACGT\n>x\nACGT\n')
    with pytest.raises(ValueError, match='duplicate header'):
        load(dup)


def test_duplicate_headers_fail_validation_exit_2(tmp_path, capsys):
    dup = write(tmp_path, 'd.fa', '>x\nACGT\n>x\nTGCA\n')
    ok = write(tmp_path, 'r.fa', '>x\nACGT\n')
    assert main([str(dup), str(ok)]) == 2
    assert 'duplicate header' in capsys.readouterr().err


def test_sequence_before_header_raises(tmp_path):
    bad = write(tmp_path, 'b.fa', 'ACGT\n>x\nACGT\n')
    with pytest.raises(ValueError, match='before any header'):
        load(bad)


# --- header id extraction ------------------------------------------------------------

def test_repbase_tab_delimited_header_uses_first_token(tmp_path):
    fa = write(tmp_path, 'r.fa', '>ALU\tSINE1/7SL\tPrimates\nACGT\n')
    assert set(load(fa)) == {'ALU'}


def test_blank_lines_ignored(tmp_path):
    a = write(tmp_path, 'a.fa', '>x\n\nACGT\n\n')
    b = write(tmp_path, 'b.fa', '>x\nACGT\n')
    assert load(a) == load(b)


# --- comparison report ---------------------------------------------------------------

def test_report_counts_all_five_categories(tmp_path):
    new = write(tmp_path, 'n.fa', '>same\nACGT\n>diff\nACGT\n>extra\nACGT\n')
    ref = write(tmp_path, 'r.fa', '>same\nACGT\n>diff\nTTTT\n>missing\nACGT\n')
    r = compare(load(new), load(ref))
    assert r['new_total'] == 3 and r['ref_total'] == 3
    assert r['shared'] == 2
    assert r['missing'] == ['MISSING']
    assert r['extra'] == ['EXTRA']
    assert r['matching'] == 1
    assert r['mismatching'] == ['DIFF']


def test_same_header_different_sequence_is_caught(tmp_path, capsys):
    """The MER5A case: identical header, different sequence must FAIL."""
    new = write(tmp_path, 'n.fa', '>MER5A\nACGTACGT\n')
    ref = write(tmp_path, 'r.fa', '>MER5A\nACGTACGA\n')
    assert main([str(new), str(ref)]) == 1
    out = capsys.readouterr().out
    assert 'mismatching=1' in out and 'MER5A' in out and 'FAIL' in out


def test_iupac_difference_does_not_count_as_mismatch(tmp_path, capsys):
    """RepBase .ref vs index differ only by IUPAC->N; canonicalization absorbs it."""
    new = write(tmp_path, 'n.fa', '>MER5A\nACGTNNGT\n')
    ref = write(tmp_path, 'r.fa', '>MER5A\nACGTRYGT\n')
    assert main([str(new), str(ref)]) == 0
    assert 'matching=1' in capsys.readouterr().out


# --- threshold enforcement -----------------------------------------------------------

def test_threshold_enforced_on_match_fraction(tmp_path, capsys):
    body = ''.join(f'>f{i}\nACGT\n' for i in range(99)) + '>f99\nACGT\n'
    new = write(tmp_path, 'n.fa', body)
    ref = write(tmp_path, 'r.fa', body.replace('>f99\nACGT\n', '>f99\nTTTT\n'))
    assert main([str(new), str(ref), '--min-match', '0.99']) == 0    # 99/100
    assert main([str(new), str(ref), '--min-match', '1.0']) == 1
    assert 'match_fraction=0.99000' in capsys.readouterr().out


def test_missing_headers_fail_even_when_shared_all_match(tmp_path):
    new = write(tmp_path, 'n.fa', '>a\nACGT\n')
    ref = write(tmp_path, 'r.fa', '>a\nACGT\n>b\nACGT\n')
    assert main([str(new), str(ref)]) == 1                       # default --max-missing 0
    assert main([str(new), str(ref), '--max-missing', '1']) == 0


def test_extra_headers_alone_do_not_fail(tmp_path):
    new = write(tmp_path, 'n.fa', '>a\nACGT\n>b\nACGT\n')
    ref = write(tmp_path, 'r.fa', '>a\nACGT\n')
    assert main([str(new), str(ref)]) == 0
