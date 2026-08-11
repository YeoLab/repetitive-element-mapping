"""Generate MASTER_FILELIST TSV (and identical .list) from Gencode, RepBase, rmsk, tRNA, miRNA and rRNA.

The file is a CONCATENATION OF SOURCE LISTS, not a table. Column 5 records which
list a row came from, and the same name may legitimately appear twice with
different annotations (61 do in the hg38 reference, e.g. DNA1_MAM as TcMar in
the RepBase block and TcMar-Tc1 in the rmsk block). Do not deduplicate globally
and do not sort.

Block order, matching the hg38 reference:

    Gencode + rRNA   ENST | ENSG | gene_name | FAMILY | genelists.FAMILY
    RepBase families NAME | FAM  | FAM       | FAM    | CLASS
    tRNA             name | anticodon | anticodon | tRNA | <gtrnadb file>
    SimpleRepeat     AT_SimpleRepeat | Simple_repeat x4
    miRNA            MI0000001 | miRNA x3 | <mirbase name>   (+ '-proximal' twin)
    rmsk leftovers   NAME | FAM | FAM | FAM | CLASS, and (XXX)N simple repeats
"""

import argparse
import collections
import re
import sys
from pathlib import Path

import repbase

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import (
    open_maybe_gz, parse_gtf_attributes, setup_logger,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_bowtie2_index import (  # noqa: E402
    RRNA_SUBUNITS, read_gtrnadb_fasta, read_rrna_genbank, read_rrna_subunit_args,
    simple_repeat_fasta,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gtf', required=True,
                   help='Gencode GTF, for gene_id / gene_name / gene_type.')
    p.add_argument('--parsed-ucsc', required=True,
                   help='parsed_ucsc_tableformat, for transcript coordinates.')
    p.add_argument('--repeatmasker', required=True,
                   help='Assembly RepeatMasker GTF. Supplies the repeat NAME set '
                        'and its order only; class and family come from '
                        '--rmsk-class-family.')
    p.add_argument('--rmsk-class-family', required=True,
                   help='TSV of repName/repClass/repFamily, e.g. '
                        'refdata/hg38.rmsk-class-family.tsv.gz. Regenerate with '
                        '`zcat rmsk.txt.gz | cut -f11,12,13 | sort -u`. The '
                        'assembly RepeatMasker GTF in examples/inputs is stripped '
                        'to gene_id and cannot supply these.')
    p.add_argument('--rmsk-smallrna-bed',
                   help='BED4 of rmsk small-RNA loci (repClass in srpRNA, scRNA, '
                        'snRNA, tRNA, rRNA, RNA), e.g. '
                        'refdata/hg38.rmsk-smallrna.bed.gz. Drives the strongest '
                        'Gencode family rule -- on hg38 it resolves 3,932 of '
                        '5,261 rows. Without it, families come from gene_name '
                        'and --family-override only.')
    p.add_argument('--repbase-species-fasta', required=True,
                   help='RepBase 18.05 species FASTA, the same input T-05 takes. '
                        'Selection is repbase.select_families, so the RepBase '
                        'block matches the index exactly without the filelist '
                        'having to wait for an index to exist -- the build order '
                        'is filelist, then index.')
    p.add_argument('--repbase-drop-exact',
                   help='Override the RepBase drop-exact list. Mouse needs '
                        'refdata/repbase18.05-drop-exact.mm.txt; the default is '
                        'the human list.')
    p.add_argument('--repbase-class-family',
                   help='Curated family -> repFamily/repClass table, '
                        'refdata/hg38.repbase-class-family.tsv. Consulted after '
                        'the rmsk table. Covers the 354 hg38 families with no '
                        'genomic instances in the modern rmsk track, and is the '
                        'cross-species source for mouse -- pass the HUMAN file '
                        'to a mouse build (-475.42).')
    p.add_argument('--family-order', required=True,
                   help='refdata/gencode-family-order.txt.')
    p.add_argument('--chrom-allowlist', required=True)
    p.add_argument('--gtrnadb-fasta')
    p.add_argument('--gff3')
    p.add_argument('--rrna-genbank', action='append', default=[])
    p.add_argument('--rrna-subunit', action='append', default=[], metavar='LABEL=PATH')
    p.add_argument('--family-override', action='append', default=[],
                   help='TSV of transcript_id<TAB>family, applied before the '
                        'rmsk and gene_name rules. Repeatable. Intended for the '
                        'residue neither rule resolves (-475.41).')
    p.add_argument('--output', required=True)
    return p.parse_args()


# ── Gencode family assignment ────────────────────────────────────────────
#
# Measured on hg38: rmsk overlap resolves 3,932/5,261, gene_name another 831,
# and 473 are resolved by neither (clone-style names like AC020634.1 whose only
# Gencode annotation is gene_type "snoRNA"). Those need --family-override.

RMSK_REPNAME_TO_FAMILY = {
    'U1': 'RNU1', 'U2': 'RNU2', 'U4': 'RNU4', 'U6': 'RNU6', 'U7': 'RNU7',
    'U11': 'RNU11', 'U12': 'RNU12', 'U4atac': 'RNU4ATAC', 'U6atac': 'RNU6ATAC',
    '7SK': 'RN7SK', '7SLRNA': 'RN7SL', '5S': 'RNA5S',
    'HY1': 'YRNA', 'HY3': 'YRNA', 'HY4': 'YRNA', 'HY5': 'YRNA',
    'U3': 'SNORD', 'U8': 'SNORD', 'U13_': 'SNORD', 'U14': 'SNORD',
    'ACA64': 'SNORA',
}

# U5 subfamilies and U17 are the two repNames that map to more than one family;
# both are resolved by gene_name instead.
RMSK_REPNAME_TO_FAMILY_UPPER = {
    k.upper(): v for k, v in RMSK_REPNAME_TO_FAMILY.items()
}

RMSK_AMBIGUOUS_REPNAMES = frozenset({'U5', 'U17'})

# Used only to report residue: a transcript of one of these types that resolves
# to no family is a real gap (-475.41). Everything else in the GTF is meant to
# be absent, so silence is correct there.
SMALL_RNA_GENE_TYPES = frozenset({
    'snRNA', 'snoRNA', 'misc_RNA', 'scaRNA', 'scRNA', 'sRNA', 'rRNA',
    'rRNA_pseudogene', 'Mt_tRNA', 'Mt_rRNA', 'vault_RNA', 'vaultRNA', 'ribozyme',
})

GENE_NAME_FAMILY_RULES = [
    (re.compile(r'^RNVU1-'), 'RNU1'),
    (re.compile(r'^RNU6ATAC'), 'RNU6ATAC'),
    (re.compile(r'^RNU4ATAC'), 'RNU4ATAC'),
    (re.compile(r'^(RNU5[A-Z])'), None),          # RNU5A..RNU5F keep the letter
    (re.compile(r'^(RNU\d+)'), None),             # RNU1-3 -> RNU1
    (re.compile(r'^(RN7SK|RN7SL)'), None),
    (re.compile(r'^(SNORD|SNORA|SCARNA)'), None),
    (re.compile(r'^(VTRNA\d+)'), None),
    (re.compile(r'^RNY'), 'YRNA'),
    (re.compile(r'^Y_RNA'), 'YRNA'),
    (re.compile(r'^(RNA5-8S|RNA5S|RNA18S|RNA28S|RNA45S)'), None),
    (re.compile(r'^MT-T'), 'MTTRNA'),
    (re.compile(r'^MT-RNR1'), 'MTRNR1'),
    (re.compile(r'^MT-RNR2'), 'MTRNR2'),
    (re.compile(r'^Mt-t', re.I), 'MTTRNA'),
]


def family_from_gene_name(gene_name):
    """Return the curated family for a gene_name, or None."""
    for pattern, fixed in GENE_NAME_FAMILY_RULES:
        m = pattern.match(gene_name)
        if not m:
            continue
        if fixed:
            return fixed
        return m.group(1).upper() if m.groups() else m.group(0).upper()
    return None


def read_family_overrides(paths):
    """Read transcript_id -> family TSVs."""
    out = {}
    for path in paths:
        for line in Path(path).read_text().splitlines():
            if not line.strip() or line.startswith('#'):
                continue
            tid, family = line.split('\t')[:2]
            out[tid.strip()] = family.strip()
    return out


def read_family_order(path):
    return [
        line.strip() for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith('#')
    ]


def read_chrom_allowlist(path):
    return {
        line.strip() for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith('#')
    }


def read_gtf_lookup(gtf_path):
    """transcript_id -> (gene_id, gene_name, gene_type)."""
    out = {}
    with open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'transcript':
                continue
            a = parse_gtf_attributes(parts[8])
            tid = a.get('transcript_id', '')
            if tid:
                out[tid] = (a.get('gene_id', ''), a.get('gene_name', ''),
                            a.get('gene_type', a.get('gene_biotype', '')))
    return out


def read_parsed_ucsc_loci(path):
    """transcript_id -> (chrom, txStart, txEnd), in file order."""
    out = collections.OrderedDict()
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 6 or p[1] in out:
                continue
            out[p[1]] = (p[2], int(p[4]), int(p[5]))
    return out


def read_rmsk_class_family(path):
    """repName (uppercased) -> (repClass, repFamily), first occurrence wins.

    Trailing '?' marks a provisional RepeatMasker call ('DNA?'); the reference
    carries the settled label, so it is stripped.

    UCSC qualifies some repNames with a subfamily after a slash -- alpha
    satellite is 'ALR/Alpha' -- while the index and the filelist use the bare
    'ALR'. The part before the slash is indexed as well, which is what recovers
    Satellite/centr for the ALR* family and 350-odd others.
    """
    out = {}
    aliases = {}
    with open_maybe_gz(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 3:
                continue
            value = (p[1].rstrip('?'), p[2].rstrip('?'))
            name = p[0].upper()
            out.setdefault(name, value)
            if '/' in name:
                aliases.setdefault(name.split('/', 1)[0], value)
    for name, value in aliases.items():
        out.setdefault(name, value)
    return out


# RepBase tags the species in the family name. Only these are species tags:
# '_LTR', '_DNA', '_II', '_MAM' and friends are class tokens inside the name
# (ERVB3_1-LTR_MM), and stripping those would corrupt the family.
SPECIES_TAGS = ('_MM', '_HS')


def resolve_repeat(name, class_family, curated_class_family):
    """Return (repFamily, repClass, how) for a repeat name.

    Resolution order, most authoritative first:

      1. the curated small-RNA family map -- 5S is RNA5S everywhere in the
         reference, never rRNA/rRNA
      2. the assembly rmsk table, including UCSC's slash-qualified aliases
      3. the curated RepBase table (-475.42), which covers the 354 hg38
         families that have no genomic instances in the modern rmsk track, and
         doubles as the cross-species source: 170 mouse families share an exact
         family name with a human one
      4. the same two tables again after stripping a species tag, which
         recovers 27 further mouse families (RLTR1_MM -> RLTR1 -> ERV1/LTR)

    rmsk is consulted before the curated table so that current RepeatMasker
    calls win and post-2020 reclassification stays visible.
    """
    upper = name.upper()
    curated_small_rna = RMSK_REPNAME_TO_FAMILY_UPPER.get(upper)
    if curated_small_rna:
        return curated_small_rna, curated_small_rna, 'small_rna'

    for candidate, how in repeat_name_aliases(upper):
        if candidate in class_family:
            rep_class, rep_family = class_family[candidate]
            return rep_family, rep_class, 'rmsk' + how
        if candidate in curated_class_family:
            return (*curated_class_family[candidate], 'curated' + how)

    return name, name, 'unresolved'


def repeat_name_aliases(upper):
    """Yield (candidate name, suffix describing how it was derived).

    The exact name is tried first, so an alias never beats a direct hit.

    Two systematic naming differences are covered. RepeatMasker writes the
    internal segment of an LTR element as NAME_I-int where RepBase writes
    NAME_I -- that recovers the highest-footprint mouse residue, ERVs carrying
    1-2 Mb each. And RepBase tags the species in the name, so RLTR1_MM is the
    same family the mouse rmsk track calls RLTR1.
    """
    def with_int(base):
        yield base, ''
        if base.endswith('_I'):
            yield base + '-INT', '_int'

    yield from with_int(upper)
    for tag in SPECIES_TAGS:
        if upper.endswith(tag):
            for candidate, how in with_int(upper[:-len(tag)]):
                yield candidate, '_stem' + how


def repeat_row(name, class_family, curated_class_family=None):
    """A repeat row: NAME | FAM | FAM | FAM | CLASS."""
    fam, cls, _ = resolve_repeat(name, class_family, curated_class_family or {})
    return (name, fam, fam, fam, cls)


def read_curated_class_family(path):
    """family -> (repFamily, repClass) from refdata/<asm>.repbase-class-family.tsv."""
    out = {}
    if not path:
        return out
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#') or line.startswith('family\t'):
            continue
        p = line.split('\t')
        if len(p) >= 3:
            out[p[0].upper()] = (p[1], p[2])
    return out


def read_rmsk_small_rna_loci(path):
    """BED4 of small-RNA loci as {chrom: [(start, end, repName), ...]} sorted."""
    by_chrom = collections.defaultdict(list)
    with open_maybe_gz(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                continue
            by_chrom[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def overlapping_repname(loci, chrom, start, end, min_frac=0.5):
    """repName of the small-RNA feature covering >= min_frac of [start, end)."""
    best, best_ov = None, 0
    span = max(1, end - start)
    for f_start, f_end, name in loci.get(chrom, ()):
        if f_start >= end:
            break
        ov = min(end, f_end) - max(start, f_start)
        if ov > best_ov:
            best, best_ov = name, ov
    return best if best_ov / span >= min_frac else None


def read_repeatmasker_names(path):
    """(repeat_names, simple_names) as UPPERCASE, in first-occurrence order."""
    repeats, simples = [], []
    seen_r, seen_s = set(), set()
    with open_maybe_gz(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'exon':
                continue
            gid = parse_gtf_attributes(p[8]).get('gene_id', '')
            if not gid:
                continue
            up = gid.upper()
            if re.match(r'^\(.+\)N$', up):
                if up not in seen_s:
                    seen_s.add(up)
                    simples.append(up)
            elif up not in seen_r:
                seen_r.add(up)
                repeats.append(up)
    return repeats, simples


def read_mirna_gff3(path):
    """(id, name) for miRNA_primary_transcript features, in file order."""
    out, seen = [], set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'miRNA_primary_transcript':
                continue
            mid = re.search(r'ID=([^;]+)', p[8])
            name = re.search(r'Name=([^;]+)', p[8])
            if mid and name and mid.group(1) not in seen:
                seen.add(mid.group(1))
                out.append((mid.group(1), name.group(1)))
    return out


# ── block builders ───────────────────────────────────────────────────────

def build_gencode_rows(args, log):
    """Gencode and rRNA rows, grouped into per-family blocks."""
    gtf = read_gtf_lookup(args.gtf)
    loci = read_parsed_ucsc_loci(args.parsed_ucsc)
    allowed = read_chrom_allowlist(args.chrom_allowlist)
    overrides = read_family_overrides(args.family_override)
    small_rna = (read_rmsk_small_rna_loci(args.rmsk_smallrna_bed)
                 if args.rmsk_smallrna_bed else {})

    by_family = collections.defaultdict(list)
    stats = collections.Counter()
    for tid, (chrom, start, end) in loci.items():
        if chrom not in allowed or tid not in gtf:
            continue
        gene_id, gene_name, gene_type = gtf[tid]
        family = overrides.get(tid)
        if family:
            stats['override'] += 1
        else:
            rep = overlapping_repname(small_rna, chrom, start, end)
            if rep and rep not in RMSK_AMBIGUOUS_REPNAMES:
                family = RMSK_REPNAME_TO_FAMILY.get(rep)
            if family:
                stats['rmsk'] += 1
            else:
                family = family_from_gene_name(gene_name)
                if family:
                    stats['gene_name'] += 1
        if not family:
            # Having no family is the selection filter, not an error: most of
            # the GTF is protein-coding and simply does not belong here. Only a
            # small-RNA transcript that resolves to nothing is real residue.
            if gene_type in SMALL_RNA_GENE_TYPES:
                stats['unresolved_small_rna'] += 1
            continue
        by_family[family].append((tid, gene_id, gene_name, family,
                                  f'genelists.{family}'))

    log.info(f'  Gencode family source: {dict(stats)}')
    order = read_family_order(args.family_order)
    rows = []
    for family in order:
        rows.extend(by_family.pop(family, []))
    for family in sorted(by_family):                 # families not in the file
        log.warning(f'  family {family} not in --family-order; appended')
        rows.extend(by_family[family])
    return rows, stats


RRNA_GENE_RE = re.compile(r'/gene="RNA45S(N\d+)"')


def rrna_copy_suffix(path):
    """The 'N5' of /gene="RNA45SN5", identifying which rDNA copy this is.

    col3 in the reference is RNA18SN5, not a positional index, so it has to be
    read from the record rather than enumerated.
    """
    m = RRNA_GENE_RE.search(Path(path).read_text())
    return m.group(1) if m else ''


def build_rrna_rows(args):
    """rRNA rows keyed like the reference: RNA18SN5 / RNA18S / genelists.RNA18S."""
    if not args.rrna_genbank:
        return []
    subunits = read_rrna_subunit_args(args.rrna_subunit)
    per_file = [(read_rrna_genbank(p, subunits), rrna_copy_suffix(p))
                for p in args.rrna_genbank]
    rows = []
    for subunit in RRNA_SUBUNITS + ('45S',):
        for records, suffix in per_file:
            for name in records:
                if name.endswith('-' + subunit):
                    family = f'RNA{subunit}'
                    rows.append((name, name, family + suffix, family,
                                 f'genelists.{family}'))
    return rows


def build_repbase_rows(species_fasta, class_family, log, drop_exact_path=None,
                       curated_class_family=None):
    """The 1,224 RepBase families: NAME | FAM | FAM | FAM | CLASS."""
    drop_exact = (repbase.load_drop_exact(drop_exact_path)
                  if drop_exact_path else None)
    families, unparsed = repbase.select_families(species_fasta, drop_exact)
    if unparsed:
        log.error(f'{len(unparsed)} RepBase records have no derivable family '
                  f'name: {unparsed[:10]}')
        sys.exit(1)
    curated_class_family = curated_class_family or {}
    rows, how = [], collections.Counter()
    unresolved = []
    for name, _, _ in families:
        name = name.upper()
        fam, cls, source = resolve_repeat(name, class_family, curated_class_family)
        how[source] += 1
        if source == 'unresolved':
            unresolved.append(name)
        rows.append((name, fam, fam, fam, cls))
    log.info(f'  RepBase class/family source: {dict(how)}')
    if unresolved:
        log.warning(f'  {len(unresolved)} RepBase families resolved by nothing; '
                    f'class/family fall back to the name itself, so each counts '
                    f'as its own family (e.g. {unresolved[:5]}) -- see -475.42')
    return rows


def build_trna_rows(gtrnadb_path):
    """Each locus followed immediately by its _withgenomeflank twin.

    col5 names the list the rows came from. The reference writes
    'hg38-tRNAs.fa.list.wgenome_flank.flankNs', describing the flank-doubled,
    N-padded list the tRNA block actually is.
    """
    source = Path(gtrnadb_path).name + '.list.wgenome_flank.flankNs'
    rows = []
    for name, *_ in read_gtrnadb_fasta(gtrnadb_path):
        anticodon = name.rsplit('-', 2)[0]
        rows.append((name, anticodon, anticodon, 'tRNA', source))
        rows.append((name + '_withgenomeflank', anticodon + '_withgenomeflank',
                     anticodon + '_withgenomeflank', 'tRNA', source))
    return rows


def build_simple_repeat_rows():
    """The same 501 synthesized k-mers the index carries."""
    return [
        (line[1:], 'Simple_repeat', 'Simple_repeat', 'Simple_repeat',
         'Simple_repeat')
        for line in simple_repeat_fasta().splitlines() if line.startswith('>')
    ]


def build_mirna_rows(gff3_path):
    rows = []
    for mid, name in read_mirna_gff3(gff3_path):
        rows.append((mid, 'miRNA', 'miRNA', 'miRNA', name))
        rows.append((mid + '-proximal', 'miRNA-proximal', 'miRNA-proximal',
                     'miRNA-proximal', name + '-proximal'))
    return rows


def build_rmsk_rows(repeat_names, simple_names, class_family, already,
                    curated_class_family=None):
    """Everything in the assembly track the earlier blocks did not cover."""
    rows = []
    for name in repeat_names:
        if name in already:
            continue
        rows.append(repeat_row(name, class_family, curated_class_family))
    for name in simple_names:
        rows.append((name, 'Simple_repeat', 'Simple_repeat', 'Simple_repeat',
                     'Simple_repeat'))
    return rows


def main():
    args = parse_args()
    log = setup_logger('generate_master_filelist')

    if not args.rmsk_smallrna_bed:
        log.warning('--rmsk-smallrna-bed not provided; Gencode families will '
                    'come from gene_name and --family-override only')

    log.info('Building Gencode rows...')
    gencode_rows, stats = build_gencode_rows(args, log)
    log.info(f'  {len(gencode_rows)} Gencode rows')

    log.info('Building rRNA rows...')
    rrna_rows = build_rrna_rows(args)
    log.info(f'  {len(rrna_rows)} rRNA rows')

    # rRNA sits inside the Gencode region, in --family-order position.
    order = read_family_order(args.family_order)
    combined = gencode_rows + rrna_rows
    rank = {f: i for i, f in enumerate(order)}
    blocks = collections.OrderedDict()
    for row in combined:
        blocks.setdefault(row[3], []).append(row)
    head_rows = []
    for family in sorted(blocks, key=lambda f: rank.get(f, len(rank))):
        head_rows.extend(blocks[family])

    curated_class_family = read_curated_class_family(args.repbase_class_family)
    if curated_class_family:
        log.info(f'Curated RepBase class/family entries: '
                 f'{len(curated_class_family)}')

    log.info('Reading rmsk class/family table...')
    class_family = read_rmsk_class_family(args.rmsk_class_family)
    log.info(f'  {len(class_family)} repeat names')

    log.info('Building RepBase rows...')
    repbase_rows = build_repbase_rows(args.repbase_species_fasta, class_family, log,
                                      args.repbase_drop_exact, curated_class_family)
    log.info(f'  {len(repbase_rows)} RepBase rows')

    trna_rows = build_trna_rows(args.gtrnadb_fasta) if args.gtrnadb_fasta else []
    if args.gtrnadb_fasta:
        log.info(f'  {len(trna_rows)} tRNA rows (plain + genome-flanked)')
    else:
        log.warning('--gtrnadb-fasta not provided; tRNA rows omitted')

    simple_rows = build_simple_repeat_rows()
    log.info(f'  {len(simple_rows)} synthesized simple-repeat rows')

    mirna_rows = build_mirna_rows(args.gff3) if args.gff3 else []
    if args.gff3:
        log.info(f'  {len(mirna_rows)} miRNA rows')
    else:
        log.warning('--gff3 not provided; miRNA rows omitted')

    log.info('Reading assembly RepeatMasker track...')
    repeat_names, simple_names = read_repeatmasker_names(args.repeatmasker)
    log.info(f'  {len(repeat_names)} repeat names, {len(simple_names)} (X)N names')
    rmsk_rows = build_rmsk_rows(repeat_names, simple_names, class_family,
                                {r[0] for r in repbase_rows}, curated_class_family)
    log.info(f'  {len(rmsk_rows)} rmsk leftover rows')

    rows = (head_rows + repbase_rows + trna_rows + simple_rows + mirna_rows
            + rmsk_rows)
    log.info(f'Total rows: {len(rows)}')

    out_path = Path(args.output)
    text = ''.join('\t'.join(r) + '\n' for r in rows)
    out_path.write_text(text)
    out_path.with_suffix('.list').write_text(text)
    log.info(f'Wrote {out_path} and {out_path.with_suffix(".list")}')

    if stats["unresolved_small_rna"]:
        log.warning(f'{stats["unresolved_small_rna"]} small-RNA transcripts had no family '
                    f'from rmsk, gene_name or --family-override and were '
                    f'omitted (see -475.41)')


if __name__ == '__main__':
    main()
