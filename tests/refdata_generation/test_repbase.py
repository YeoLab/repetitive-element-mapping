"""RepBase family selection and NAME extraction (issue -7ee / T-05).

The load-bearing claim: NAME extraction reproduces all 1,224 hg38 repeat-family
names with zero curated exceptions. hg38 is the only ground truth there is, so
if these break, the mouse build has nothing to stand on.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin" / "python" / "refdata_generation"))
import repbase  # noqa: E402

HUMAN = Path(
    "/tscc/projects/ps-yeolab4/genomes/RepBase18.05.fasta/species_specific/"
    "homo_sapiens_repbase_fixed_v2.fasta"
)
MOUSE = Path(
    "/tscc/projects/ps-yeolab4/genomes/RepBase18.05.fasta/species_specific/"
    "mus_musculus_repbase_u1_fixed_v2.fastq"
)
INDEX = REPO / (
    "examples/inputs/hg38/bowtie2_index/"
    "MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa"
)

needs_repbase = pytest.mark.skipif(
    not HUMAN.exists(), reason="RepBase 18.05 not available"
)
needs_index = pytest.mark.skipif(
    not INDEX.exists(), reason="hg38 reference index not available"
)


# --- taxon stripping ----------------------------------------------------------

@pytest.mark.parametrize("header,expected", [
    ("ALUY_SINE1/7SL_Homo_sapiens", "ALUY_SINE1/7SL"),   # binomial taxon
    ("MER5A_hAT_Mammalia", "MER5A_hAT"),                  # single-token taxon
    ("ALRa__SAT_Primates", "ALRa__SAT"),                  # empty token survives
])
def test_strip_taxon(header, expected):
    assert repbase.strip_taxon(header) == expected


# --- NAME extraction ----------------------------------------------------------

def test_longest_class_suffix_wins():
    """CR1_HS_CR1_Homo_sapiens is NAME=CR1_HS CLASS=CR1.

    A left-to-right class match returns 'CR1' -- the wrong family.
    """
    classes = {"CR1", "hAT"}
    assert repbase.family_name("CR1_HS_CR1_Homo_sapiens", classes) == "CR1_HS"


def test_bare_header_is_its_own_name():
    """The mouse file stores U-snRNAs as bare headers with no class or taxon."""
    assert repbase.family_name("U1", {"snRNA"}) == "U1"
    assert repbase.family_name("UHG", {"snRNA"}) == "UHG"


def test_trailing_underscore_name_preserved():
    """ALRa__SAT_Primates has NAME 'ALRa_' -- an empty second token."""
    assert repbase.family_name("ALRa__SAT_Primates", {"SAT"}) == "ALRa_"


# --- selection rule -----------------------------------------------------------

@pytest.mark.parametrize("header,kept", [
    ("(A)n_Simple_Repeat_Eukaryota", False),
    ("tRNA-Ala-AGC_tRNA_Homo_sapiens", False),
    ("MamSINE1_tRNA_Mammalia", True),      # tRNA-derived SINE: a real family
    ("ALUY_SINE1/7SL_Primates", True),
    ("U1_snRNA_Homo_sapiens", False),      # in DROP_EXACT: Gencode supplies it
    ("U3_snRNA_Homo_sapiens", True),       # NOT dropped: a blanket snRNA drop is wrong
])
def test_selection_rule(header, kept):
    assert repbase.keep_record(header, repbase.load_drop_exact()) is kept


# --- the claim that everything else rests on ----------------------------------

@needs_repbase
@needs_index
def test_human_reproduces_reference_index_exactly():
    families, unparsed = repbase.select_families(HUMAN)
    index_headers = {
        line[1:].strip().upper() for line in INDEX.open() if line.startswith(">")
    }
    names = {name for name, _, _ in families}

    assert unparsed == [], f"unparsed headers: {unparsed[:5]}"
    assert len(families) == 1224, f"selection gave {len(families)}, expected 1224"
    missing = names - index_headers
    assert not missing, f"{len(missing)} names absent from the index: {sorted(missing)[:5]}"


@needs_repbase
def test_mouse_duplicate_family_is_fatal():
    """U7/U14/U8 exist twice in the mouse file -- bare and class-formatted, with
    byte-identical sequences. Emitting both would put one family in the index
    twice, so this must raise rather than silently keep one."""
    with pytest.raises(ValueError, match="duplicate family"):
        repbase.select_families(MOUSE)
