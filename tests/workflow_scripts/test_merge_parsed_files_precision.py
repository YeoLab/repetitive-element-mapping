"""RPR must be written at full precision (issue -475.23).

merge_parsed_files.py formatted the reads-per-read column with ':.5f'. The Perl
original writes a bare double. On the se_full SE dataset that truncation drove
39 of 182 elements' Input_clip_rpr to exactly 0.0, which makes Fold_enrichment
'inf' -- and where the IP side also truncated to zero, 0/0 left the cell empty.
86 cells were destroyed and 636 float columns were wrong by up to 37%.

calculate_fold_change_from_parsed_files.py reads clip_rpr straight out of this
column, so any rounding here lands directly in the pipeline's primary output.
"""

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "workflow" / "scripts" / "merge_parsed_files.py"

# A count small enough that count/usable underflows 5 decimal places.
# 12 / 17_691_900 = 6.78e-07, which ':.5f' renders as "0.00000".
RARE_COUNT = 12
USABLE = 17_691_900


def write_parsed(path, total_lines, usable=USABLE):
    """Minimal per-prefix .parsed file in the format merge_parsed_files reads."""
    body = [
        f"#READINFO\tAllReads\t{usable}",
        f"#READINFO\tUsableReads\t{usable}\t1.0",
        f"#READINFO\tGenomicReads\t0\t0.0",
        f"#READINFO\tRepFamilyReads\t{usable}\t1.0",
    ]
    body += [f"TOTAL\t{el}\t{n}\t{n / usable}" for el, n in total_lines]
    path.write_text("\n".join(body) + "\n")
    return path


def run_merge(tmp_path, inputs):
    out = tmp_path / "merged.parsed"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(out), *(str(i) for i in inputs)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return out


def rpr_for(merged, element):
    for line in merged.read_text().splitlines():
        parts = line.split("\t")
        if parts[0] == "TOTAL" and parts[1] == element:
            return parts[3]
    raise AssertionError(f"{element} not found in {merged}")


def test_rare_element_rpr_does_not_truncate_to_zero(tmp_path):
    """The exact failure: a rare element's RPR must not become 0.00000."""
    src = write_parsed(tmp_path / "AA.parsed", [("antisense_Crypton", RARE_COUNT)])
    merged = run_merge(tmp_path, [src])

    rpr = rpr_for(merged, "antisense_Crypton")
    assert float(rpr) != 0.0, (
        f"RPR truncated to zero ({rpr!r}); a zero denominator makes "
        "Fold_enrichment inf downstream"
    )
    assert float(rpr) == RARE_COUNT / USABLE


def test_rpr_round_trips_exactly(tmp_path):
    """Full float64 precision, not a fixed number of decimals."""
    counts = [("RNA28S", 5_399_580), ("unique_distintron", 2_443_909),
              ("rare", 1)]
    src = write_parsed(tmp_path / "AA.parsed", counts)
    merged = run_merge(tmp_path, [src])

    for element, n in counts:
        assert float(rpr_for(merged, element)) == n / USABLE, (
            f"{element} RPR lost precision"
        )


def test_readinfo_genomic_fraction_not_truncated(tmp_path):
    """The #READINFO fraction columns carried the same ':.5f' defect.

    Uses a genomic count whose fraction does not terminate within 5 decimals,
    so a truncated render is unequal to the exact quotient.
    """
    genomic = 8_214_042
    src = tmp_path / "AA.parsed"
    src.write_text(
        f"#READINFO\tAllReads\t{USABLE}\n"
        f"#READINFO\tUsableReads\t{USABLE}\t1.0\n"
        f"#READINFO\tGenomicReads\t{genomic}\t{genomic / USABLE}\n"
        f"#READINFO\tRepFamilyReads\t{USABLE - genomic}\t{(USABLE - genomic) / USABLE}\n"
        f"TOTAL\tRNA28S\t5399580\t{5399580 / USABLE}\n"
    )
    merged = run_merge(tmp_path, [src])

    fractions = {
        parts[1]: parts[3]
        for parts in (l.split("\t") for l in merged.read_text().splitlines())
        if parts[0] == "#READINFO" and len(parts) >= 4
    }
    assert float(fractions["GenomicReads"]) == genomic / USABLE, (
        f"GenomicReads fraction truncated: {fractions['GenomicReads']!r} "
        f"!= {genomic / USABLE!r}"
    )
