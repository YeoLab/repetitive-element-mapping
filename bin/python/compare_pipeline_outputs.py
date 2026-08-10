#!/usr/bin/env python3
"""Compare two pipeline .nopipes.tsv/.withpipes.tsv outputs.

Read counts are integers produced by deterministic counting and must match
exactly. The derived float columns (clip_rpr, Fold_enrichment,
Information_content) are compared with a relative tolerance, because the CWL and
Snakemake paths reach them through different arithmetic and differ at the
~1e-15 level.

Row ORDER is not compared. Elements with equal read counts are ordered by hash
iteration in the Perl implementation, so tie order is not reproducible across
implementations and carries no meaning.

Exit 0 if equivalent, 1 if not.

    compare_pipeline_outputs.py NEW.nopipes.tsv REFERENCE.nopipes.tsv
"""

import argparse
import csv
import sys

COUNT_COLS = ["IP_read_num", "Input_read_num"]
FLOAT_COLS = [
    "IP_clip_rpr",
    "Input_clip_rpr",
    "Fold_enrichment",
    "Information_content",
]


def load(path):
    with open(path, newline="") as fh:
        rows = {r["element"]: r for r in csv.DictReader(fh, delimiter="\t")}
    if not rows:
        sys.exit(f"{path}: no data rows")
    return rows


def to_float(s):
    """Parse a cell, tolerating the empty/inf values a degraded run can emit.

    A zero denominator yields 'inf'; 0/0 yields an empty cell. Both are real
    output states worth reporting rather than crashing on.
    """
    s = s.strip() if s is not None else ""
    if s == "":
        return None
    return float(s)


def rel_diff(a, b):
    """Relative difference, falling back to absolute when b is 0."""
    if b == 0:
        return abs(a)
    return abs(a - b) / abs(b)


def compare(new, ref, tol, show):
    problems = []

    only_new = sorted(set(new) - set(ref))
    only_ref = sorted(set(ref) - set(new))
    if only_new:
        problems.append(f"{len(only_new)} elements only in new: {only_new[:show]}")
    if only_ref:
        problems.append(f"{len(only_ref)} elements only in reference: {only_ref[:show]}")

    shared = sorted(set(new) & set(ref))
    print(f"elements: new={len(new)} reference={len(ref)} shared={len(shared)}")

    count_mismatch = []
    for el in shared:
        for col in COUNT_COLS:
            # Input_read_num is written as a float ("905719.0") by the
            # fold-change script even though it is a count.
            a, b = to_float(new[el][col]), to_float(ref[el][col])
            if a != b:
                count_mismatch.append(f"{el}.{col}: {a} vs {b}")
    if count_mismatch:
        problems.append(
            f"{len(count_mismatch)} read-count mismatches (must be exact): "
            + "; ".join(count_mismatch[:show])
        )
    else:
        print(f"read counts: all {len(shared)} elements match exactly")

    worst = (0.0, None)
    float_mismatch = []
    degenerate = []
    for el in shared:
        for col in FLOAT_COLS:
            a, b = to_float(new[el][col]), to_float(ref[el][col])
            if a is None or b is None or a != a or b != b or a in (
                float("inf"), float("-inf")
            ) or b in (float("inf"), float("-inf")):
                # Empty / inf / nan: not comparable as a number. Report the
                # side that degraded rather than silently skipping.
                degenerate.append(f"{el}.{col}: {new[el][col]!r} vs {ref[el][col]!r}")
                continue
            d = rel_diff(a, b)
            if d > worst[0]:
                worst = (d, f"{el}.{col} ({a} vs {b})")
            if d > tol:
                float_mismatch.append(f"{el}.{col}: {a} vs {b} (rel {d:.2e})")
    if degenerate:
        problems.append(
            f"{len(degenerate)} empty/inf/nan float cells (a zero denominator "
            f"destroys fold enrichment): " + "; ".join(degenerate[:show])
        )
    if float_mismatch:
        problems.append(
            f"{len(float_mismatch)} float columns exceed tol={tol:.1e}: "
            + "; ".join(float_mismatch[:show])
        )
    print(f"derived floats: max relative deviation {worst[0]:.3e} at {worst[1]}")

    return problems


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("new")
    p.add_argument("reference")
    p.add_argument("--tol", type=float, default=1e-9,
                   help="relative tolerance for derived float columns "
                        "(default 1e-9; CWL-vs-Snakemake noise is ~1e-15)")
    p.add_argument("--show", type=int, default=5,
                   help="how many examples to print per problem class")
    args = p.parse_args()

    problems = compare(load(args.new), load(args.reference), args.tol, args.show)

    if problems:
        print("\nFAIL")
        for pr in problems:
            print(f"  - {pr}")
        return 1
    print("\nPASS — equivalent within tolerance (row order ignored)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
