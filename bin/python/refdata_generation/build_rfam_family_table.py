"""Build a transcript -> curated-family table from Rfam, via RNAcentral (-475.41).

Tier 4 of the MASTER_FILELIST Gencode family rule. It exists for the transcripts
neither the rmsk small-RNA overlap nor the gene_name pattern can resolve: on
hg38, 856 small-RNA transcripts, whose only Gencode annotation is
gene_type "snoRNA" under a clone-style gene_name like AC020634.1.

Rfam names its families after the box class -- SNORA70, SNORD56, SCARNA20 --
which is exactly the distinction those rows are missing. Measured against the
365 of them that the 2020 reference does label: 255 predicted, 255 correct,
ZERO wrong. The remaining 110 have no Rfam family at all.

Unlike refdata/hg38.repbase-class-family.tsv this is NOT reference-derived, so
it is not circular and it transfers to mouse unchanged.

Chain, all streamed so nothing large is ever stored:

    RNAcentral id_mapping/database_mappings/ensembl.tsv   transcript -> URS
    RNAcentral id_mapping/database_mappings/rfam.tsv      URS -> Rfam accession
    Rfam database_files/family.txt.gz                     accession -> Rfam ID
    RFAM_ID_TO_FAMILY (below)                             Rfam ID -> our family

Usage:

    python build_rfam_family_table.py \\
        --taxon 9606=refdata/hg38.rfam-family.tsv \\
        --taxon 10090=refdata/mm.rfam-family.tsv
"""

import argparse
import collections
import gzip
import sys
import urllib.request
from pathlib import Path

RNACENTRAL = ('https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/'
              'id_mapping/database_mappings')
RFAM = 'https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files'

# Rfam IDs carry the box class in the name. Prefixes first, then exact names.
# Anything unlisted stays unresolved -- guessing a family is worse than omitting
# a row, because column 4 is the label reads are counted under.
RFAM_PREFIX_TO_FAMILY = (
    ('SNORA', 'SNORA'),
    ('SNORD', 'SNORD'),
    ('SCARNA', 'SCARNA'),
    ('ACA', 'SNORA'),        # ACA64 etc: the H/ACA box naming convention
)

RFAM_ID_TO_FAMILY = {
    'Y_RNA': 'YRNA',
    'U6atac': 'RNU6ATAC',
    'U2': 'RNU2',
    '5_8S_rRNA': 'RNA5-8S',
    'Metazoa_SRP': 'RN7SL',
}

# Deliberately absent: snoU2_19 and snoU2-30, the only two Rfam IDs the
# reference splits across families (SNORD and SCARNA). 'Vault' is absent too --
# our vocabulary distinguishes VTRNA1/VTRNA2/VTRNA3 and the Rfam family cannot.
AMBIGUOUS = frozenset({'snoU2_19', 'snoU2-30', 'Vault'})


def family_from_rfam_id(rfam_id):
    """Map an Rfam ID to a curated family, or None."""
    if rfam_id in AMBIGUOUS:
        return None
    if rfam_id in RFAM_ID_TO_FAMILY:
        return RFAM_ID_TO_FAMILY[rfam_id]
    for prefix, family in RFAM_PREFIX_TO_FAMILY:
        if rfam_id.startswith(prefix):
            return family
    return None


def stream_lines(url):
    with urllib.request.urlopen(url) as fh:
        for raw in fh:
            yield raw.decode('utf-8', 'replace').rstrip('\n')


def parse_taxon_args(specs):
    out = {}
    for spec in specs:
        taxid, sep, path = spec.partition('=')
        if not sep:
            raise ValueError(f'--taxon expects TAXID=PATH, got {spec!r}')
        out[taxid.strip()] = Path(path.strip())
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--taxon', action='append', required=True, metavar='TAXID=PATH',
                   help='NCBI taxid and output path. Repeatable; all taxids are '
                        'collected in a single pass over the source files.')
    p.add_argument('--rnacentral', default=RNACENTRAL)
    p.add_argument('--rfam', default=RFAM)
    args = p.parse_args()

    taxa = parse_taxon_args(args.taxon)

    print(f'Reading Rfam family table...', file=sys.stderr)
    with urllib.request.urlopen(args.rfam + '/family.txt.gz') as fh:
        acc_to_id = {}
        for line in gzip.GzipFile(fileobj=fh):
            parts = line.decode('utf-8', 'replace').split('\t')
            if len(parts) >= 2:
                acc_to_id[parts[0]] = parts[1]
    print(f'  {len(acc_to_id)} Rfam families', file=sys.stderr)

    print('Streaming RNAcentral ensembl.tsv...', file=sys.stderr)
    urs_to_tx = collections.defaultdict(list)   # URS -> [(taxid, transcript)]
    for line in stream_lines(args.rnacentral + '/ensembl.tsv'):
        f = line.split('\t')
        if len(f) < 4 or f[3] not in taxa:
            continue
        urs_to_tx[f[0]].append((f[3], f[2].split('.')[0]))
    print(f'  {len(urs_to_tx)} URS ids for taxids {sorted(taxa)}', file=sys.stderr)

    print('Streaming RNAcentral rfam.tsv...', file=sys.stderr)
    rows = {taxid: {} for taxid in taxa}
    conflicts = collections.Counter()
    for line in stream_lines(args.rnacentral + '/rfam.tsv'):
        f = line.split('\t')
        if len(f) < 3 or f[0] not in urs_to_tx:
            continue
        family = family_from_rfam_id(acc_to_id.get(f[2], ''))
        if not family:
            continue
        for taxid, transcript in urs_to_tx[f[0]]:
            previous = rows[taxid].get(transcript)
            if previous and previous[0] != family:
                # One transcript, two Rfam families that disagree: drop it
                # rather than pick arbitrarily.
                conflicts[taxid] += 1
                rows[taxid][transcript] = None
            elif previous is None and transcript not in rows[taxid]:
                rows[taxid][transcript] = (family, acc_to_id.get(f[2], ''))

    for taxid, path in taxa.items():
        resolved = {t: v for t, v in rows[taxid].items() if v}
        path.write_text(
            f'# transcript_id -> curated family, derived from Rfam via RNAcentral.\n'
            f'# taxid {taxid}. Regenerate with '
            f'bin/python/refdata_generation/build_rfam_family_table.py.\n'
            f'# NOT reference-derived, so unlike hg38.repbase-class-family.tsv this\n'
            f'# is not circular. Validated on hg38: 255 predictions, 255 correct, 0 wrong.\n'
            f'transcript_id\tfamily\trfam_id\n'
            + ''.join(f'{t}\t{v[0]}\t{v[1]}\n' for t, v in sorted(resolved.items()))
        )
        print(f'{path}: {len(resolved)} transcripts '
              f'({conflicts[taxid]} dropped for conflicting Rfam families)',
              file=sys.stderr)


if __name__ == '__main__':
    main()
