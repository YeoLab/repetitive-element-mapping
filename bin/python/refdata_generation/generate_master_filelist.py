"""Generate MASTER_FILELIST TSV (and identical .list) from Gencode, RepeatMasker, tRNA, miRNA, custom FASTA."""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import open_maybe_gz, parse_gtf_attributes, setup_logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gtf', required=True, help='Gencode GTF (for gene_id, gene_name, transcript_type lookup)')
    p.add_argument('--parsed-ucsc', required=True)
    p.add_argument('--repeatmasker', required=True)
    p.add_argument('--simplerepeats', help='Not used for FASTA; kept for CLI compatibility')
    p.add_argument('--trna')
    p.add_argument('--gff3')
    p.add_argument('--custom-fasta', action='append', default=[])
    p.add_argument('--output', required=True)
    return p.parse_args()


def _is_simple_repeat_gene_id(gid):
    return bool(re.match(r'^\(.+\)n$', gid))


def _simple_repeat_to_fasta_name(gid):
    inner = gid[1:-2]
    return inner.upper() + '_SimpleRepeat'


def read_gtf_lookup(gtf_path):
    """Return dicts: transcript_id → (gene_id, gene_name, transcript_type)."""
    tid_to_info = {}
    with open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'transcript':
                continue
            attrs = parse_gtf_attributes(parts[8])
            tid = attrs.get('transcript_id', '')
            gid = attrs.get('gene_id', '')
            gname = attrs.get('gene_name', '')
            ttype = attrs.get('transcript_type', attrs.get('transcript_biotype', ''))
            if tid:
                tid_to_info[tid] = (gid, gname, ttype)
    return tid_to_info


def read_parsed_ucsc_tids(path):
    """Return list of (gene_id_col0, transcript_id_col1) in file order."""
    pairs = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            gid_col, tid = parts[0], parts[1]
            if tid not in seen:
                seen.add(tid)
                pairs.append((gid_col, tid))
    return pairs


def read_repeatmasker_for_filelist(path):
    """Return (non_simple_rows, simple_rows).
    non_simple_rows: list of UPPERCASE gene_id (deduplicated, first occurrence)
    simple_rows: list of PATTERN_SimpleRepeat name (deduplicated)
    """
    seen_nonsimple = {}  # uppercase → True
    seen_simple = {}     # pattern name → True
    nonsimple_names = []
    simple_names = []

    with open_maybe_gz(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'exon':
                continue
            attrs = parse_gtf_attributes(parts[8])
            gid = attrs.get('gene_id', '')
            if not gid:
                continue
            if _is_simple_repeat_gene_id(gid):
                name = _simple_repeat_to_fasta_name(gid)
                if name not in seen_simple:
                    seen_simple[name] = True
                    simple_names.append(name)
            else:
                up = gid.upper()
                if up not in seen_nonsimple:
                    seen_nonsimple[up] = True
                    nonsimple_names.append(up)

    return nonsimple_names, simple_names


def read_trna_gene_ids(path):
    """Return list of unique gene_ids from tRNA tsv (deduplicated, first occurrence)."""
    seen = {}
    ids = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'exon':
                continue
            attrs = parse_gtf_attributes(parts[8])
            gid = attrs.get('gene_id', '')
            if gid and gid not in seen:
                seen[gid] = True
                ids.append(gid)
    return ids


def read_mirna_gff3_for_filelist(path):
    """Return list of (id, name) for miRNA_primary_transcript features (deduplicated)."""
    seen = {}
    entries = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'miRNA_primary_transcript':
                continue
            attrs_str = parts[8]
            id_m = re.search(r'ID=([^;]+)', attrs_str)
            name_m = re.search(r'Name=([^;]+)', attrs_str)
            if not id_m or not name_m:
                continue
            entry_id = id_m.group(1)
            entry_name = name_m.group(1)
            if entry_id not in seen:
                seen[entry_id] = True
                entries.append((entry_id, entry_name))
    return entries


def read_custom_fasta_headers(paths):
    """Return list of FASTA sequence IDs (first word after '>') from all custom FASTA files."""
    headers = []
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith('>'):
                    headers.append(line[1:].split()[0])
    return headers


DEFAULT_TRANSCRIPT_TYPES = {
    'snRNA', 'misc_RNA', 'snoRNA', 'rRNA_pseudogene', 'scaRNA',
    'rRNA', 'Mt_tRNA', 'Mt_rRNA', 'vaultRNA',
}


def main():
    args = parse_args()
    log = setup_logger('generate_master_filelist')

    # 1. GTF lookup: transcript_id → (gene_id, gene_name, transcript_type)
    log.info('Reading GTF for transcript metadata...')
    tid_to_info = read_gtf_lookup(args.gtf)
    log.info(f'  {len(tid_to_info)} transcripts in GTF')

    # 2. parsed_ucsc → ordered list of (gene_id_col0, transcript_id)
    log.info('Reading parsed_ucsc...')
    ucsc_pairs = read_parsed_ucsc_tids(args.parsed_ucsc)
    log.info(f'  {len(ucsc_pairs)} transcripts in parsed_ucsc')

    rows = []

    # Group 1: Gencode (filtered by transcript_type)
    log.info('Building Gencode rows...')
    gc_count = 0
    seen_tids = set()
    for gid_col0, tid in ucsc_pairs:
        if tid in seen_tids:
            continue
        info = tid_to_info.get(tid)
        if info is None:
            continue
        gene_id, gene_name, ttype = info
        if ttype not in DEFAULT_TRANSCRIPT_TYPES:
            continue
        seen_tids.add(tid)
        rows.append((tid, gene_id, gene_name, ttype, f'genelists.{ttype}'))
        gc_count += 1
    log.info(f'  {gc_count} Gencode rows')

    # Group 2 & 3: RepeatMasker non-simple and simple repeats
    log.info('Reading RepeatMasker...')
    rm_names, simple_names = read_repeatmasker_for_filelist(args.repeatmasker)
    log.info(f'  {len(rm_names)} non-simple, {len(simple_names)} simple repeat entries')

    for name in rm_names:
        rows.append((name, name, name, name, name))

    for name in simple_names:
        rows.append((name, 'Simple_repeat', 'Simple_repeat', 'Simple_repeat', 'Simple_repeat'))

    # Group 4: tRNA
    if args.trna:
        log.info('Reading tRNA...')
        trna_ids = read_trna_gene_ids(args.trna)
        log.info(f'  {len(trna_ids)} tRNA entries')
        for gid in trna_ids:
            rows.append((gid, gid, gid, 'tRNA', 'genelists.tRNA'))
    else:
        log.warning('--trna not provided; tRNA entries omitted')

    # Group 5: miRNA + miRNA-proximal
    if args.gff3:
        log.info('Reading miRNA gff3...')
        mirna_entries = read_mirna_gff3_for_filelist(args.gff3)
        log.info(f'  {len(mirna_entries)} miRNA entries')
        for entry_id, entry_name in mirna_entries:
            rows.append((entry_id, 'miRNA', 'miRNA', 'miRNA', entry_name))
            rows.append((entry_id + '-proximal', 'miRNA-proximal', 'miRNA-proximal',
                         'miRNA-proximal', entry_name + '-proximal'))
    else:
        log.warning('--gff3 not provided; miRNA entries omitted')

    # Group 6: Custom FASTA (rRNA)
    if args.custom_fasta:
        log.info(f'Reading {len(args.custom_fasta)} custom FASTA file(s)...')
        headers = read_custom_fasta_headers(args.custom_fasta)
        log.info(f'  {len(headers)} custom FASTA headers')
        for h in headers:
            rows.append((h, h, h, 'rRNA', 'genelists.rRNA'))

    # Deduplicate by col1 (sequence_id), keep first occurrence
    seen_ids = set()
    deduped = []
    for row in rows:
        if row[0] not in seen_ids:
            seen_ids.add(row[0])
            deduped.append(row)

    log.info(f'Total rows (after dedup): {len(deduped)}')

    # Write TSV
    out_path = Path(args.output)
    log.info(f'Writing {out_path}...')
    with open(out_path, 'w') as out:
        for row in deduped:
            out.write('\t'.join(row) + '\n')

    # Write identical .list sibling
    list_path = out_path.with_suffix('.list')
    log.info(f'Writing {list_path}...')
    with open(list_path, 'w') as out:
        for row in deduped:
            out.write('\t'.join(row) + '\n')

    log.info('Done.')


if __name__ == '__main__':
    main()
