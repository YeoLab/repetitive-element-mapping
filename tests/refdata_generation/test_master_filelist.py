"""MASTER_FILELIST generation (T-09, -gz2).

The file is a concatenation of source lists, not a table. Two premises in the
original T-09 scope were wrong and both are pinned here by tests:

  * repeat class/family are the UCSC rmsk repClass/repFamily, NOT the middle
    field of the RepBase header -- L2B_CR1_Eutheria is family L2, class LINE.
  * rows must NOT be deduplicated by column 1; the hg38 reference carries 61
    names twice with different annotations.
"""

import gzip
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin" / "python" / "refdata_generation"))

REFERENCE = REPO / (
    "examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA."
    "enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv"
)
FAMILY_ORDER = REPO / "refdata/gencode-family-order.txt"

needs_reference = pytest.mark.skipif(
    not REFERENCE.exists(), reason="hg38 reference MASTER_FILELIST not available"
)


# ── Gencode family from gene_name ────────────────────────────────────────

@pytest.mark.parametrize("gene_name, expected", [
    ("RNU1-3", "RNU1"),
    ("RNVU1-22", "RNU1"),          # RNVU1 is a U1 variant, not its own family
    ("RNU5A-1", "RNU5A"),          # U5 keeps its subfamily letter
    ("RNU5F-1", "RNU5F"),
    ("RNU6-2", "RNU6"),
    ("RNU6ATAC", "RNU6ATAC"),      # must win over the plain RNU6 rule
    ("RNU4ATAC", "RNU4ATAC"),
    ("RN7SL1", "RN7SL"),
    ("SNORD3A", "SNORD"),
    ("SNORA73B", "SNORA"),
    ("SCARNA5", "SCARNA"),
    ("VTRNA1-1", "VTRNA1"),
    ("RNY1", "YRNA"),
    ("MT-TA", "MTTRNA"),
    ("MT-RNR1", "MTRNR1"),
    ("MT-RNR2", "MTRNR2"),
    ("AC020634.1", None),          # clone-style: needs --family-override
    ("ACTB", None),
])
def test_family_from_gene_name(gene_name, expected):
    from generate_master_filelist import family_from_gene_name
    assert family_from_gene_name(gene_name) == expected


def test_rnu6atac_is_not_swallowed_by_the_rnu6_rule():
    """Rule order matters: ^(RNU\\d+) would otherwise claim RNU6ATAC."""
    from generate_master_filelist import family_from_gene_name
    assert family_from_gene_name("RNU6ATAC") != "RNU6"


# ── repeat rows ──────────────────────────────────────────────────────────

def _class_family_file(tmp_path, rows):
    p = tmp_path / "cf.tsv.gz"
    with gzip.open(p, "wt") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return p


def test_repeat_row_uses_rmsk_family_and_class(tmp_path):
    """L2B is family L2 / class LINE, though RepBase calls its class CR1."""
    from generate_master_filelist import read_rmsk_class_family, repeat_row
    cf = read_rmsk_class_family(
        _class_family_file(tmp_path, [("L2b", "LINE", "L2")]))
    assert repeat_row("L2B", cf) == ("L2B", "L2", "L2", "L2", "LINE")


def test_provisional_question_mark_is_stripped(tmp_path):
    from generate_master_filelist import read_rmsk_class_family, repeat_row
    cf = read_rmsk_class_family(
        _class_family_file(tmp_path, [("MER105", "DNA?", "hAT-Charlie?")]))
    assert repeat_row("MER105", cf)[4] == "DNA"
    assert repeat_row("MER105", cf)[1] == "hAT-Charlie"


def test_slash_qualified_repname_is_reachable_by_its_bare_prefix(tmp_path):
    """UCSC writes alpha satellite as ALR/Alpha; the filelist uses ALR."""
    from generate_master_filelist import read_rmsk_class_family, repeat_row
    cf = read_rmsk_class_family(
        _class_family_file(tmp_path, [("ALR/Alpha", "Satellite", "centr")]))
    assert repeat_row("ALR", cf) == ("ALR", "centr", "centr", "centr",
                                     "Satellite")


def test_exact_name_wins_over_a_slash_alias(tmp_path):
    from generate_master_filelist import read_rmsk_class_family, repeat_row
    cf = read_rmsk_class_family(_class_family_file(tmp_path, [
        ("ALR/Alpha", "Satellite", "centr"),
        ("ALR", "Satellite", "exact"),
    ]))
    assert repeat_row("ALR", cf)[1] == "exact"


def test_small_rna_repnames_carry_the_curated_family(tmp_path):
    """5S is RNA5S everywhere in the reference, never rRNA/rRNA."""
    from generate_master_filelist import read_rmsk_class_family, repeat_row
    cf = read_rmsk_class_family(
        _class_family_file(tmp_path, [("5S", "rRNA", "rRNA")]))
    assert repeat_row("5S", cf) == ("5S", "RNA5S", "RNA5S", "RNA5S", "RNA5S")


def test_unknown_repeat_falls_back_to_its_own_name(tmp_path):
    from generate_master_filelist import read_rmsk_class_family, repeat_row
    cf = read_rmsk_class_family(_class_family_file(tmp_path, []))
    assert repeat_row("ALR1", cf) == ("ALR1",) * 5


# ── tRNA and rRNA blocks ─────────────────────────────────────────────────

def test_trna_rows_interleave_each_locus_with_its_flank_twin(tmp_path):
    from generate_master_filelist import build_trna_rows
    fa = tmp_path / "hg38-tRNAs.fa"
    fa.write_text(
        ">Homo_sapiens_tRNA-Ala-AGC-1-1 x chr6:100-171 (-)\nACGT\n"
        ">Homo_sapiens_tRNA-Gly-CCC-2-1 x chr1:10-81 (+)\nACGT\n"
    )
    rows = build_trna_rows(fa)
    assert [r[0] for r in rows] == [
        "tRNA-Ala-AGC-1-1", "tRNA-Ala-AGC-1-1_withgenomeflank",
        "tRNA-Gly-CCC-2-1", "tRNA-Gly-CCC-2-1_withgenomeflank",
    ]
    assert rows[0][1] == "tRNA-Ala-AGC"          # anticodon level, copy dropped
    assert rows[0][3] == "tRNA"
    assert rows[0][4] == "hg38-tRNAs.fa.list.wgenome_flank.flankNs"


def test_rrna_copy_suffix_comes_from_the_gene_qualifier(tmp_path):
    """col3 is RNA18SN5, identifying the rDNA copy -- not a positional index."""
    from generate_master_filelist import rrna_copy_suffix
    gb = tmp_path / "r.gb"
    gb.write_text('VERSION     NR_046235.3\n     gene            1..10\n'
                  '                     /gene="RNA45SN5"\n')
    assert rrna_copy_suffix(gb) == "N5"


# ── structural invariants of the reference ───────────────────────────────

@needs_reference
def test_reference_repeats_column_one_so_dedup_would_be_wrong():
    """61 names appear twice with DIFFERENT annotations. Never dedup by col1."""
    import collections
    col1 = collections.Counter(l.split("\t")[0] for l in REFERENCE.open())
    repeated = {k for k, n in col1.items() if n > 1}
    assert len(repeated) == 61
    rows = collections.defaultdict(set)
    for line in REFERENCE.open():
        p = line.rstrip("\n").split("\t")
        if p[0] in repeated:
            rows[p[0]].add(tuple(p))
    assert any(len(v) > 1 for v in rows.values()), (
        "duplicates must differ in annotation, else dedup would be harmless"
    )


@needs_reference
def test_family_order_file_matches_the_reference_block_order():
    from generate_master_filelist import read_family_order
    seen = []
    for line in REFERENCE.open():
        p = line.rstrip("\n").split("\t")
        if not (p[0].startswith("ENST") or p[0].startswith("NR_")):
            break
        if not seen or seen[-1] != p[3]:
            seen.append(p[3])
    assert read_family_order(FAMILY_ORDER) == seen
