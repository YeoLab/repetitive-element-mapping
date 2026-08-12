"""Cross-artifact consistency checks for a reference-data set (issue -475.13 / M-8).

Every hg38 acceptance criterion is "reproduce the reference file". Mouse has no
reference, so correctness cannot be established that way. These six checks
replace reference-diff with internal consistency between the three artifacts of
one assembly -- MASTER_FILELIST, bowtie2 index FASTA, UniqueGenomicElements BED
-- and they are calibrated so that they PASS on the hg38 reference set and FAIL
on each of the three known-bad artifacts this project has produced:

    T-05 (index built per genomic instance)   -> check C, then A
    T-07 (BED over-selected from 5 sources)   -> check B, then E
    M-5  (RepBase block duplicates Gencode)   -> check D

Exit status: 0 if every non-skipped check passes, 1 if any fails, 2 on
malformed input.

Two of the six cannot be stated the way the issue words them, and the wording
is wrong rather than the artifacts:

  A. "every MASTER_FILELIST id resolves to an index header AND VICE VERSA."
     Only the index->filelist direction is an invariant. 18,748 of hg38's
     26,354 filelist ids (71%) have no index sequence by design -- the rmsk
     leftover and SimpleRepeat blocks contribute NAMES, which the dedup perl
     needs in order to type a peak, not sequences to align against. Requiring
     the reverse direction would fail on the reference set.

  D. "no family appears in both the repeat and Gencode portions."
     The hg38 reference itself has two, SNORD and YRNA, from the RepBase
     U3/U8/U13/U14 snoRNA records it deliberately keeps. So the check is that
     the overlap is a SUBSET of the reference's own overlap, i.e. mouse may not
     invent an overlap human does not have. Mouse currently sits at {SNORD}.
"""

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import setup_logger  # noqa: E402


# The RepBase-vs-Gencode family overlap the hg38 REFERENCE set carries. Anything
# in an assembly's overlap beyond this is a defect; anything less is fine.
REFERENCE_FAMILY_OVERLAP = frozenset({'SNORD', 'YRNA'})

# Column 5 of the tRNA block, e.g. 'hg38-tRNAs.fa.list.wgenome_flank.flankNs'.
# Used to find where the RepBase block ends. Column 4 cannot do it: mm10's
# MamSINE1 row carries family 'tRNA' (mm39's carries 'tRNA-RTE'), taken from
# each assembly's own rmsk table, and would cut the block 684 rows early.
TRNA_SOURCE_SUFFIX = '.list.wgenome_flank.flankNs'

# RepeatMasker rows are ~99.8% of every BED. A set where they are not the bulk
# has either lost them or gained a source that does not belong (the T-07
# simple-repeat over-generation put them at 78%).
MIN_RMSK_FRACTION = 0.99


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--assembly', required=True)
    p.add_argument('--master-filelist', required=True)
    p.add_argument('--index-fasta', required=True)
    p.add_argument('--unique-genomic-elements', required=True)
    p.add_argument('--repbase-provenance',
                   help='<index prefix>.repbase_provenance.tsv. Check F is '
                        'skipped when absent -- the hg38 reference index is a '
                        'downloaded 2020 artifact and has no sidecar.')
    p.add_argument('--show', type=int, default=10,
                   help='offending items to list per failure (default 10)')
    return p.parse_args()


def read_filelist(path):
    """Return (ids, gencode_families, repbase_families, repbase_names).

    ids: every column-1 id, uppercased, pipe-joined ids split. This is the key
    space read_in_filelists builds convert_enst2type from, so it is what an
    index header or a BED name has to hit.

    The file is a concatenation of source-list blocks in a fixed order (see
    docs/MASTER_FILELIST-decisions.md), so the RepBase block is the contiguous
    run between the last genelists row and the first tRNA row.
    """
    rows = [line.rstrip('\n').split('\t') for line in open(path) if line.strip()]
    ids = set()
    for r in rows:
        if r[0]:
            ids.update(r[0].upper().split('|'))

    gencode_idx = [i for i, r in enumerate(rows)
                   if len(r) >= 5 and r[4].startswith('genelists.')]
    trna_idx = [i for i, r in enumerate(rows)
                if len(r) >= 5 and r[4].endswith(TRNA_SOURCE_SUFFIX)]
    if not gencode_idx or not trna_idx:
        raise ValueError(f'{path}: cannot locate the Gencode or tRNA block '
                         f'(genelists rows={len(gencode_idx)}, tRNA rows={len(trna_idx)})')

    gencode_families = {rows[i][3] for i in gencode_idx}
    repbase_block = rows[max(gencode_idx) + 1:min(trna_idx)]
    # Column 4 is the FAMILY (SINE, SNORD); column 1 is the family NAME the
    # index headers and the provenance sidecar are keyed by (B1_MM, ALINE).
    repbase_families = {r[3] for r in repbase_block if len(r) >= 4}
    repbase_names = {r[0].upper() for r in repbase_block if r[0]}
    return ids, gencode_families, repbase_families, repbase_names


def read_index_headers(path):
    """FASTA header ids, uppercased, in file order."""
    return [line[1:].strip().upper() for line in open(path) if line.startswith('>')]


def read_bed_sources(path, filelist_ids):
    """Return (total, unresolved_rows, per-source counts, distinct miRNA ids).

    Source is inferred from the row's own shape, the same way the generator
    writes it: Gencode carries '-' in column 5, tRNA carries a tRNA-<aa>-<anti>
    family, miRNA is named MI<digits>, the trf simple-repeat track is named
    trf*, and everything else is RepeatMasker.
    """
    counts = collections.Counter()
    mirna_ids = set()
    total = unresolved = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            total += 1
            if f[3].upper() not in filelist_ids:
                unresolved += 1
            if f[4] == '-':
                counts['gencode'] += 1
            elif f[4].startswith('tRNA-'):
                counts['trna'] += 1
            elif f[3].startswith('MI0'):
                counts['mirna'] += 1
                mirna_ids.add(f[3].split('-proximal')[0])
            elif f[3].startswith('trf'):
                counts['simplerepeat'] += 1
            else:
                counts['rmsk'] += 1
    return total, unresolved, counts, mirna_ids


class Report:
    """Collects per-check outcomes so every check runs before anything exits."""

    def __init__(self, log, show):
        self.log = log
        self.show = show
        self.failed = 0
        self.skipped = 0

    def _items(self, items):
        shown = sorted(items)[:self.show]
        extra = len(items) - len(shown)
        return ', '.join(map(str, shown)) + (f' ... (+{extra} more)' if extra else '')

    def check(self, key, name, ok, detail, offenders=()):
        if ok:
            self.log.info(f'  {key} PASS  {name}: {detail}')
        else:
            self.failed += 1
            self.log.error(f'  {key} FAIL  {name}: {detail}')
            if offenders:
                self.log.error(f'         offending: {self._items(offenders)}')

    def skip(self, key, name, why):
        self.skipped += 1
        self.log.warning(f'  {key} SKIP  {name}: {why}')


def main():
    args = parse_args()
    log = setup_logger('validate_refdata_set')
    log.info(f'Validating {args.assembly}')

    try:
        ids, gencode_families, repbase_families, repbase_names = read_filelist(args.master_filelist)
        headers = read_index_headers(args.index_fasta)
        bed_total, bed_unresolved, bed_counts, mirna_ids = read_bed_sources(
            args.unique_genomic_elements, ids)
    except (OSError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2

    log.info(f'  filelist {len(ids)} ids ({len(repbase_names)} RepBase families, '
             f'{len(gencode_families)} Gencode families) | index {len(headers)} headers | '
             f'BED {bed_total} rows')

    r = Report(log, args.show)

    # A. index headers must all be typeable. An index sequence whose header is
    # absent from the filelist aligns reads that read_in_filelists cannot type.
    orphans = [h for h in headers if h not in ids]
    r.check('A', 'index headers resolve in MASTER_FILELIST',
            not orphans, f'{len(headers) - len(orphans)}/{len(headers)}', orphans)

    # B. same contract for the BED: read_peakfi uppercases column 4 and looks it
    # up in convert_enst2type, so an unresolved name is an untyped peak.
    r.check('B', 'BED names resolve in MASTER_FILELIST',
            bed_unresolved == 0,
            f'{bed_total - bed_unresolved}/{bed_total} rows', ())

    # C. '::' is pybedtools getfasta's coordinate suffix -- its presence means
    # the repeat portion was built per genomic instance, which was T-05.
    coord_headers = [h for h in headers if '::' in h]
    r.check('C', "no index header carries a '::' coordinate suffix",
            not coord_headers, f'{len(coord_headers)} of {len(headers)}', coord_headers)

    # D. a family in both portions splits its reads between them.
    overlap = gencode_families & repbase_families
    extra = overlap - REFERENCE_FAMILY_OVERLAP
    r.check('D', 'RepBase/Gencode family overlap is within the hg38 reference set',
            not extra,
            f'overlap={sorted(overlap) or "none"} '
            f'(reference carries {sorted(REFERENCE_FAMILY_OVERLAP)})', extra)

    # E. BED source profile. The simple-repeat track and the RepeatMasker share
    # are the two that moved under T-07; the miRNA identity is structural --
    # the generator emits one row per entry plus two proximal twins. Requiring
    # tRNA and miRNA to be non-empty is not pedantry: every optional source is
    # a --flag the generator merely WARNS about, and mm39's superseded BED had
    # zero of both because neither input existed yet.
    rmsk_fraction = bed_counts['rmsk'] / bed_total if bed_total else 0.0
    profile = ' '.join(f'{k}={v}' for k, v in sorted(bed_counts.items()))
    failures = []
    if bed_counts['simplerepeat']:
        failures.append(f"{bed_counts['simplerepeat']} trf simple-repeat rows (must be 0)")
    if rmsk_fraction < MIN_RMSK_FRACTION:
        failures.append(f'rmsk is {rmsk_fraction:.5f} of rows (min {MIN_RMSK_FRACTION})')
    if bed_counts['mirna'] != 3 * len(mirna_ids):
        failures.append(f"miRNA rows {bed_counts['mirna']} != 3 x {len(mirna_ids)} entries")
    for source in ('trna', 'mirna', 'gencode'):
        if not bed_counts[source]:
            failures.append(f'no {source} rows at all -- was the source omitted?')
    r.check('E', 'BED source profile',
            not failures,
            f'{profile} | rmsk={rmsk_fraction:.5f} | '
            + ('; '.join(failures) if failures else 'all source rules hold'))

    # F. every repeat family the index emits must be traceable to the RepBase
    # record it came from. Guards the T-05 failure mode from the other side.
    if not args.repbase_provenance:
        r.skip('F', 'provenance sidecar covers the RepBase block',
               'no --repbase-provenance given')
    else:
        covered = set()
        with open(args.repbase_provenance) as fh:
            next(fh, None)
            for line in fh:
                if line.strip():
                    covered.add(line.split('\t')[0].upper())
        indexed_repbase = {h for h in headers if h in repbase_names}
        missing = indexed_repbase - covered
        r.check('F', 'provenance sidecar covers the RepBase block',
                not missing,
                f'{len(indexed_repbase - missing)}/{len(indexed_repbase)} indexed '
                f'RepBase families have a provenance row', missing)

    verdict = 'PASS' if r.failed == 0 else 'FAIL'
    log.info(f'{args.assembly}: {verdict} '
             f'({6 - r.failed - r.skipped} passed, {r.failed} failed, {r.skipped} skipped)')
    return 0 if r.failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
