"""Compare two FASTAs by per-header canonical sequence digest (issue -475.5).

Header-level overlap is not sufficient: two FASTAs can share a header and still hold
different sequences (MER5A differs between RepBase 24.01 and the RepeatMasker Edition
while keeping the same name). This validator canonicalizes each record and compares
SHA-256 digests, streaming one record at a time so a whole-genome FASTA never has to be
held in memory.

Canonicalization: uppercase FIRST, then map every non-ACGT letter to N. Doing it in the
other order leaves lowercase ambiguity codes (rymk) as RYMK instead of NNNN.

Exit status: 0 if the shared-header match fraction meets --min-match and there are no
missing headers beyond --max-missing; 1 otherwise; 2 on malformed input.
"""

import argparse
import hashlib
import sys

# Everything that is not A/C/G/T becomes N once uppercased: the IUPAC ambiguity codes
# RYMKSWBDHV, plus N itself and any other stray letter.
_KEEP = frozenset(b'ACGT')


def canonical_digest(chunks):
    """SHA-256 of the uppercased, non-ACGT->N sequence formed by concatenating chunks."""
    h = hashlib.sha256()
    for chunk in chunks:
        upper = chunk.upper()
        h.update(bytes(b if b in _KEEP else ord('N') for b in upper))
    return h.hexdigest()


def iter_records(path):
    """Yield (header_id, digest) per FASTA record. header_id is the first whitespace-
    or tab-delimited token after '>', uppercased (RepBase headers are NAME<TAB>CLASS)."""
    header = None
    chunks = []
    with open(path, 'rb') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b'>'):
                if header is not None:
                    yield header, canonical_digest(chunks)
                header = line[1:].split()[0].decode('utf-8', 'replace').upper() if line[1:].split() else ''
                chunks = []
            else:
                if header is None:
                    raise ValueError(f'{path}: sequence data before any header line')
                chunks.append(line)
    if header is not None:
        yield header, canonical_digest(chunks)


def load(path):
    """Return {header: digest}, raising on duplicate headers."""
    out = {}
    dupes = []
    for header, digest in iter_records(path):
        if header in out:
            dupes.append(header)
        out[header] = digest
    if dupes:
        raise ValueError(
            f'{path}: {len(dupes)} duplicate header(s), e.g. '
            + ', '.join(sorted(set(dupes))[:5])
        )
    return out


def compare(new, ref):
    """Return a report dict comparing {header: digest} maps."""
    new_keys, ref_keys = set(new), set(ref)
    shared = new_keys & ref_keys
    matching = {k for k in shared if new[k] == ref[k]}
    return {
        'new_total': len(new_keys),
        'ref_total': len(ref_keys),
        'shared': len(shared),
        'missing': sorted(ref_keys - new_keys),     # in reference, absent from new
        'extra': sorted(new_keys - ref_keys),       # in new, absent from reference
        'matching': len(matching),
        'mismatching': sorted(shared - matching),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('new_fasta')
    p.add_argument('reference_fasta')
    p.add_argument('--min-match', type=float, default=0.99,
                   help='required matching/shared fraction (default 0.99)')
    p.add_argument('--max-missing', type=int, default=0,
                   help='reference headers allowed to be absent from new (default 0)')
    p.add_argument('--show', type=int, default=20,
                   help='how many offending headers to list per category (default 20)')
    args = p.parse_args(argv)

    try:
        new = load(args.new_fasta)
        ref = load(args.reference_fasta)
    except ValueError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    r = compare(new, ref)
    frac = r['matching'] / r['shared'] if r['shared'] else 0.0

    print(f"new={r['new_total']} reference={r['ref_total']} shared={r['shared']}")
    print(f"missing={len(r['missing'])} extra={len(r['extra'])}")
    print(f"matching={r['matching']} mismatching={len(r['mismatching'])} "
          f"match_fraction={frac:.5f} (threshold {args.min_match})")
    for label in ('missing', 'extra', 'mismatching'):
        if r[label]:
            shown = r[label][:args.show]
            print(f"  {label}: " + ', '.join(shown)
                  + (f' ... (+{len(r[label]) - len(shown)} more)' if len(r[label]) > len(shown) else ''))

    ok = frac >= args.min_match and len(r['missing']) <= args.max_missing
    print('PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
