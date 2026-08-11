"""tRNA and rRNA portions of the bowtie2 index (-475.36, -475.37).

Both portions were previously believed to be blocked on unobtainable source
data. They are not:

tRNA  -- the reference was built from gtRNAdb's own FASTA (hg38-tRNAs.fa), not
         from the assembly tRNA track (hg38.trna.tsv.gz). Using the track gave
         different names and a 0/413 sequence match, which is what produced the
         "different gtRNAdb release" conclusion. The gtRNAdb FASTA reproduces
         all 432 names and, under two fixed transforms, all 864 records.

rRNA  -- all 15 NR_ records come from five RefSeq GenBank flat files. The 45S
         is the whole record; 18S and 28S are cut from the misc_feature spans
         whose /note names them, so no coordinate is ever hardcoded.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin" / "python" / "refdata_generation"))

GTRNADB = REPO / "examples/inputs/hg38/downloaded/hg38-tRNAs.fa"
RRNA_DIR = REPO / "examples/inputs/hg38/downloaded/rrna"
IDX = REPO / (
    "examples/inputs/hg38/bowtie2_index/"
    "MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa"
)

needs_trna = pytest.mark.skipif(
    not (GTRNADB.exists() and IDX.exists()),
    reason="hg38 gtRNAdb FASTA / reference index not available",
)
needs_rrna = pytest.mark.skipif(
    not (RRNA_DIR.exists() and IDX.exists()),
    reason="hg38 rRNA GenBank records / reference index not available",
)


def _reference_records():
    from generate_bowtie2_index import _iter_fasta
    return {h.split()[0]: s for h, s in _iter_fasta(IDX.read_text())}


# ── tRNA ─────────────────────────────────────────────────────────────────

def test_gtrnadb_header_yields_name_and_locus(tmp_path):
    from generate_bowtie2_index import read_gtrnadb_fasta
    f = tmp_path / "t.fa"
    f.write_text(
        ">Homo_sapiens_tRNA-Ala-AGC-1-1 (tRNAscan-SE ID: chr6.trna116) Ala "
        "(AGC) 72 bp Sc: 84.9 chr6:28795964-28796035 (-)\nacgt\nACGT\n"
    )
    assert read_gtrnadb_fasta(f) == [
        ("tRNA-Ala-AGC-1-1", "ACGTACGT", "chr6", 28795963, 28796035, "-"),
    ]


def test_gtrnadb_locus_is_normalized_when_header_lists_it_descending(tmp_path):
    """Minus-strand loci are written high-to-low in some gtRNAdb releases."""
    from generate_bowtie2_index import read_gtrnadb_fasta
    f = tmp_path / "t.fa"
    f.write_text(">Homo_sapiens_tRNA-Ala-AGC-1-1 x chr6:200-100 (-)\nAC\n")
    _, _, _, start0, end0, _ = read_gtrnadb_fasta(f)[0]
    assert (start0, end0) == (99, 200)


def test_unparseable_gtrnadb_header_is_an_error(tmp_path):
    from generate_bowtie2_index import read_gtrnadb_fasta
    f = tmp_path / "t.fa"
    f.write_text(">Homo_sapiens_tRNA-Ala-AGC-1-1 no coordinates here\nAC\n")
    with pytest.raises(ValueError, match="Unparseable"):
        read_gtrnadb_fasta(f)


@needs_trna
def test_plain_trna_records_match_reference_exactly():
    """plain = N*10 + gtRNAdb genomic sequence + CCA + N*10."""
    from generate_bowtie2_index import read_gtrnadb_fasta, TRNA_PAD, TRNA_CCA
    ref = _reference_records()
    entries = read_gtrnadb_fasta(GTRNADB)
    assert len(entries) == 432
    mismatched = [
        name for name, seq, *_ in entries
        if ref.get(name) != TRNA_PAD + seq + TRNA_CCA + TRNA_PAD
    ]
    assert mismatched == []


@needs_trna
def test_every_trna_has_a_genomeflank_twin_in_the_reference():
    from generate_bowtie2_index import read_gtrnadb_fasta, TRNA_PAD, TRNA_CCA
    ref = _reference_records()
    entries = read_gtrnadb_fasta(GTRNADB)
    assert all(name + "_withgenomeflank" in ref for name, *_ in entries)
    # The flank record carries 50 bp of genomic context per side and no CCA,
    # so it is a constant 2*50 - len('CCA') longer than its plain twin.
    deltas = {
        len(ref[name + "_withgenomeflank"]) - len(ref[name])
        for name, *_ in entries
    }
    assert deltas == {2 * 50 - len(TRNA_CCA)}


# ── rRNA ─────────────────────────────────────────────────────────────────

MINIMAL_GB = """\
LOCUS       NR_000001               12 bp    rRNA    linear   PRI 01-JAN-2020
VERSION     NR_000001.1
FEATURES             Location/Qualifiers
     rRNA            1..12
                     /gene="RNA45SN0"
     misc_feature    3..5
                     /note="18S rRNA"
     misc_feature    6..7
                     /note="5.8S rRNA"
     misc_feature    9..12
                     /note="28S rRNA"
ORIGIN
        1 aacccggtta cg
//
"""


def test_genbank_subunits_are_cut_from_annotated_spans(tmp_path):
    from generate_bowtie2_index import read_rrna_genbank
    f = tmp_path / "r.gb"
    f.write_text(MINIMAL_GB)
    assert read_rrna_genbank(f) == {
        "NR_000001.1-45S": "AACCCGGTTACG",
        "NR_000001.1-18S": "CCC",
        "NR_000001.1-28S": "TACG",
    }


def test_5_8s_is_annotated_but_deliberately_not_emitted(tmp_path):
    """The reference index carries only 18S, 28S and 45S."""
    from generate_bowtie2_index import read_rrna_genbank
    f = tmp_path / "r.gb"
    f.write_text(MINIMAL_GB)
    assert not any("5.8S" in k or "5-8S" in k for k in read_rrna_genbank(f))


def test_missing_subunit_annotation_is_an_error(tmp_path):
    from generate_bowtie2_index import read_rrna_genbank
    f = tmp_path / "r.gb"
    f.write_text(MINIMAL_GB.replace('/note="28S rRNA"', '/note="ETS"'))
    with pytest.raises(ValueError, match=r"28S"):
        read_rrna_genbank(f)


# ── rRNA: unannotated precursors (mouse) ─────────────────────────────────

def _genbank(accession, seq):
    """Render a GenBank flat file with no misc_feature, like NR_046233.2."""
    body = "\n".join(
        f"{i + 1:>9} {seq[i:i + 60]}" for i in range(0, len(seq), 60)
    )
    return (
        f"LOCUS       {accession.split('.')[0]}  {len(seq)} bp rRNA\n"
        f"VERSION     {accession}\n"
        "FEATURES             Location/Qualifiers\n"
        f"     rRNA            1..{len(seq)}\n"
        '                     /product="45S pre-ribosomal RNA"\n'
        f"ORIGIN\n{body}\n//\n"
    )


def _random_seq(n, seed):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_anchors_emit_precursor_bases_not_the_subunit_record(tmp_path):
    """The mouse 28S record and the precursor's copy differ by a 3 bp indel.

    The emitted record must be the precursor's own bases over the anchored
    span, so a whole-substring match is not required and must not be assumed.
    """
    from generate_bowtie2_index import read_rrna_genbank
    head, middle, tail = _random_seq(40, 1), _random_seq(60, 2), _random_seq(40, 3)
    subunit = head + middle + tail
    in_precursor = head + middle[:30] + middle[33:] + tail      # 3 bp deleted
    precursor = _random_seq(100, 4) + in_precursor + _random_seq(100, 5)

    gb = tmp_path / "p.gb"
    gb.write_text(_genbank("NR_000002.2", precursor))
    recs = read_rrna_genbank(gb, {"18S": subunit, "28S": subunit})

    assert recs["NR_000002.2-45S"] == precursor
    assert recs["NR_000002.2-18S"] == in_precursor
    assert len(recs["NR_000002.2-18S"]) == len(subunit) - 3


def test_unannotated_subunit_without_a_fallback_is_an_error(tmp_path):
    from generate_bowtie2_index import read_rrna_genbank
    gb = tmp_path / "p.gb"
    gb.write_text(_genbank("NR_000002.2", _random_seq(300, 6)))
    with pytest.raises(ValueError, match=r"--rrna-subunit 18S"):
        read_rrna_genbank(gb)


def test_ambiguous_anchor_is_an_error(tmp_path):
    """rDNA is internally repetitive; a k-mer hitting twice must not be guessed."""
    from generate_bowtie2_index import read_rrna_genbank
    subunit = _random_seq(100, 7)
    precursor = subunit + _random_seq(50, 8) + subunit
    gb = tmp_path / "p.gb"
    gb.write_text(_genbank("NR_000002.2", precursor))
    with pytest.raises(ValueError, match="matches the precursor 2 times"):
        read_rrna_genbank(gb, {"18S": subunit, "28S": subunit})


def test_subunit_too_short_to_anchor_is_an_error(tmp_path):
    from generate_bowtie2_index import read_rrna_genbank
    gb = tmp_path / "p.gb"
    gb.write_text(_genbank("NR_000002.2", _random_seq(300, 9)))
    with pytest.raises(ValueError, match="too short to anchor"):
        read_rrna_genbank(gb, {"18S": "ACGT" * 10, "28S": "ACGT" * 10})


@pytest.mark.parametrize("spec, match", [
    ("NR_003278.3.fasta", "LABEL=PATH"),
    ("5.8S=x.fasta", "label must be one of"),
])
def test_rrna_subunit_arg_validation(spec, match):
    from generate_bowtie2_index import read_rrna_subunit_args
    with pytest.raises(ValueError, match=match):
        read_rrna_subunit_args([spec])


def test_rrna_subunit_file_must_hold_one_record(tmp_path):
    from generate_bowtie2_index import read_rrna_subunit_args
    f = tmp_path / "two.fasta"
    f.write_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(ValueError, match="expected 1 FASTA record, found 2"):
        read_rrna_subunit_args([f"18S={f}"])


MOUSE_RRNA = REPO / "examples/inputs/mm10/downloaded/rrna"

needs_mouse_rrna = pytest.mark.skipif(
    not (MOUSE_RRNA / "NR_046233.2.gb").exists(),
    reason="mouse rRNA GenBank record not available",
)


@needs_mouse_rrna
def test_mouse_precursor_yields_three_records_cut_from_itself():
    from generate_bowtie2_index import read_rrna_genbank, read_rrna_subunit_args
    subs = read_rrna_subunit_args([
        f"18S={MOUSE_RRNA / 'NR_003278.3.fasta'}",
        f"28S={MOUSE_RRNA / 'NR_003279.1.fasta'}",
    ])
    recs = read_rrna_genbank(MOUSE_RRNA / "NR_046233.2.gb", subs)
    assert set(recs) == {
        "NR_046233.2-45S", "NR_046233.2-18S", "NR_046233.2-28S",
    }
    precursor = recs["NR_046233.2-45S"]
    assert len(precursor) == 13400
    # 18S is an exact copy of Rn18s; 28S is 3 bp shorter than Rn28s1.
    assert len(recs["NR_046233.2-18S"]) == 1870
    assert len(recs["NR_046233.2-28S"]) == 4727
    assert all(recs[k] in precursor for k in recs if not k.endswith("-45S"))


@needs_rrna
def test_rrna_records_match_reference_exactly():
    from generate_bowtie2_index import rrna_fasta, _iter_fasta
    ref = _reference_records()
    built = dict(_iter_fasta(rrna_fasta(sorted(RRNA_DIR.glob("*.gb")))))
    assert len(built) == 15
    assert {k for k in ref if k.startswith("NR_")} == set(built)
    assert [k for k, v in built.items() if ref[k] != v] == []
