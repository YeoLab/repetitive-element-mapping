"""Gencode portion is selected by MASTER_FILELIST membership, not type (-475.35).

The type filter could never reproduce the reference: the 706 transcripts it
over-selected span the same types as the reference, and the reference itself
contains lncRNA, miRNA and unprocessed_pseudogene entries. The real rule is
"in the MASTER_FILELIST and on an allowlisted chromosome" -- the 259 filelist
ids the reference omits are all on GenBank scaffolds (GL000251.2, KZ208915.1).
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin" / "python" / "refdata_generation"))

ML = REPO / (
    "examples/inputs/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA."
    "enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list"
)
PU = REPO / (
    "examples/inputs/hg38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf"
    ".parsed_ucsc_tableformat"
)
IDX = REPO / (
    "examples/inputs/hg38/bowtie2_index/"
    "MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa"
)
ALLOW = REPO / "refdata/hg38.chrom-allowlist.txt"

needs_hg38 = pytest.mark.skipif(
    not (ML.exists() and PU.exists() and IDX.exists()),
    reason="hg38 reference data not available",
)


def test_read_master_filelist_splits_pipe_joined_ids(tmp_path):
    from generate_bowtie2_index import read_master_filelist_ensts
    f = tmp_path / "ml.tsv"
    f.write_text(
        "ENST00000001.1\tENSG1\tG1\tFAM\tgenelists.FAM\n"
        "ENST00000002.1|ENST00000003.1\tENSG2\tG2\tFAM\tgenelists.FAM\n"
        "ALUY\tALUY\tALUY\tALUY\tALUY\n"          # repeat row, not a transcript
    )
    assert read_master_filelist_ensts(f) == {
        "ENST00000001.1", "ENST00000002.1", "ENST00000003.1"
    }


@needs_hg38
def test_selection_reproduces_reference_transcript_set_exactly():
    from generate_bowtie2_index import (
        read_chrom_allowlist, read_master_filelist_ensts, read_parsed_ucsc,
    )
    filelist = read_master_filelist_ensts(ML)
    allowed = read_chrom_allowlist(ALLOW)
    transcripts = read_parsed_ucsc(PU)
    selected = {
        tid for tid, v in transcripts.items()
        if tid in filelist and v[0] in allowed
    }
    reference = {
        line[1:].strip() for line in IDX.open()
        if line.startswith(">") and line[1:].startswith("ENST")
    }
    assert len(selected) == 5002
    assert selected == reference, (
        f"{len(selected ^ reference)} transcripts differ from the reference"
    )
