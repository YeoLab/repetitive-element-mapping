"""SimpleRepeat portion is synthesized, not extracted (issue -7ee / T-05).

Extracting these from RepeatMasker gave 14,180 entries against the reference's
501, because it kept every observed pattern up to 10+ nt. The reference set is
exactly the primitive canonical k-mers for k=1..6 under rotation and
reverse-complement, each tiled to 60 bp.

Being pure combinatorics, this portion is species-independent -- mouse gets the
identical 501 entries with no reference to validate against.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
GEN = REPO / "bin" / "python" / "refdata_generation" / "generate_bowtie2_index.py"
REF = REPO / (
    "examples/inputs/hg38/bowtie2_index/"
    "MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa"
)


def _load_pure_functions():
    """Load only the k-mer helpers, avoiding the module's pybedtools import."""
    src = GEN.read_text()
    start = src.index("SIMPLE_REPEAT_MAX_K")
    end = src.index("def canonicalize")
    ns = {}
    exec("from itertools import product\n" + src[start:end], ns)
    return ns


FN = _load_pure_functions()


def _reference_simple_repeats():
    recs, name, seq = {}, None, []
    for line in REF.open():
        if line.startswith(">"):
            if name and "_SimpleRepeat" in name:
                recs[name] = "".join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name and "_SimpleRepeat" in name:
        recs[name] = "".join(seq)
    return recs


def _parse(fasta):
    lines = fasta.splitlines()
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 2)}


def test_primitive_rejects_repeated_units():
    assert FN["_is_primitive"]("A")
    assert FN["_is_primitive"]("AT")
    assert not FN["_is_primitive"]("AA")      # == (A)n
    assert not FN["_is_primitive"]("ATAT")    # == (AT)n


def test_canonical_kmer_folds_rotation_and_revcomp():
    canon = FN["_canonical_kmer"]
    assert canon("AT") == canon("TA")          # rotation
    assert canon("AAC") == canon("GTT")        # reverse-complement


def test_every_sequence_is_the_pattern_tiled_to_60bp():
    for header, seq in _parse(FN["simple_repeat_fasta"]()).items():
        pattern = header.replace("_SimpleRepeat", "")
        assert len(seq) == 60, f"{header} is {len(seq)}bp"
        assert seq == pattern * (60 // len(pattern))


@pytest.mark.skipif(not REF.exists(), reason="hg38 reference index not available")
def test_reproduces_reference_simple_repeats_exactly():
    mine = _parse(FN["simple_repeat_fasta"]())
    ref = _reference_simple_repeats()
    assert len(mine) == 501
    assert set(mine) == set(ref), "SimpleRepeat name sets differ"
    assert mine == ref, "SimpleRepeat sequences differ"
