"""Generate UniqueGenomicElements BED from RepeatMasker, tRNA, miRNA, and Gencode sources."""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import open_maybe_gz, parse_gtf_attributes, setup_logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repeatmasker', required=True)
    p.add_argument('--simplerepeats')
    p.add_argument('--trna')
    p.add_argument('--gff3')
    p.add_argument('--parsed-ucsc')
    p.add_argument('--assembly', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--flank', type=int, default=500)
    return p.parse_args()


def parse_gtf_bed_rows(path, name_field='gene_id'):
    """Parse GTF-format .tsv.gz; yield (chrom, start0, end0, name, score, strand) per exon."""
    with open_maybe_gz(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
            if feature != 'exon':
                continue
            attrs = parse_gtf_attributes(parts[8])
            name = attrs.get(name_field, '')
            if not name:
                continue
            start0 = int(start_s) - 1
            end0 = int(end_s)
            try:
                score_val = str(int(float(score))) if score != '.' else '0'
            except (ValueError, OverflowError):
                score_val = '0'
            yield (chrom, start0, end0, name, score_val, strand)


def parse_gff3_mirna(path):
    """Parse miRNA_primary_transcript features from gff3. Yields (chrom, start0, end0, id, name, strand)."""
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
            if feature != 'miRNA_primary_transcript':
                continue
            attrs_str = parts[8]
            id_m = re.search(r'ID=([^;]+)', attrs_str)
            name_m = re.search(r'Name=([^;]+)', attrs_str)
            if not id_m or not name_m:
                continue
            entry_id = id_m.group(1)
            entry_name = name_m.group(1)
            start0 = int(start_s) - 1
            end0 = int(end_s)
            yield (chrom, start0, end0, entry_id, entry_name, strand)


def parse_parsed_ucsc(path):
    """Yield (chrom, txStart, txEnd, transcript_id, '-', strand) per transcript."""
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            _, tid, chrom, strand, txstart, txend = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            yield (chrom, int(txstart), int(txend), tid, '-', strand)


def main():
    args = parse_args()
    log = setup_logger('generate_unique_genomic_elements')

    rows = []

    # RepeatMasker (required) — all instances, gene_id
    log.info('Reading RepeatMasker...')
    rm_rows = list(parse_gtf_bed_rows(args.repeatmasker, name_field='gene_id'))
    log.info(f'  {len(rm_rows)} RepeatMasker entries')
    rows.extend(rm_rows)

    # Simple repeats (optional) — all instances, transcript_id (Edge Case 9)
    if args.simplerepeats:
        log.info('Reading simple repeats...')
        sr_rows = list(parse_gtf_bed_rows(args.simplerepeats, name_field='transcript_id'))
        log.info(f'  {len(sr_rows)} simple repeat entries')
        rows.extend(sr_rows)
    else:
        log.warning(f'--simplerepeats not provided; simple repeat entries omitted from {args.assembly} UniqueGenomicElements')

    # tRNA (optional) — all instances, gene_id
    if args.trna:
        log.info('Reading tRNA...')
        trna_rows = list(parse_gtf_bed_rows(args.trna, name_field='gene_id'))
        log.info(f'  {len(trna_rows)} tRNA entries')
        rows.extend(trna_rows)
    else:
        log.warning(f'--trna not provided; tRNA entries omitted from {args.assembly} UniqueGenomicElements')

    # Gencode transcripts (optional) — transcript_id, score="-", actual strand
    if args.parsed_ucsc:
        log.info('Reading Gencode from parsed_ucsc...')
        gc_rows = list(parse_parsed_ucsc(args.parsed_ucsc))
        log.info(f'  {len(gc_rows)} Gencode transcript entries')
        rows.extend(gc_rows)
    else:
        log.warning(f'--parsed-ucsc not provided; Gencode entries omitted from {args.assembly} UniqueGenomicElements')

    # miRNA (optional) — ID as name, Name as score, plus 2 proximal rows per entry
    if args.gff3:
        log.info('Reading miRNA gff3...')
        mirna_count = 0
        for chrom, start0, end0, entry_id, entry_name, strand in parse_gff3_mirna(args.gff3):
            rows.append((chrom, start0, end0, entry_id, entry_name, strand))
            rows.append((chrom, max(0, start0 - args.flank), start0,
                         entry_id + '-proximal', entry_name + '-proximal', strand))
            rows.append((chrom, end0, end0 + args.flank,
                         entry_id + '-proximal', entry_name + '-proximal', strand))
            mirna_count += 1
        log.info(f'  {mirna_count} miRNA entries ({mirna_count * 2} proximal rows added)')
    else:
        log.warning(f'--gff3 not provided; miRNA entries omitted from {args.assembly} UniqueGenomicElements')

    # Sort by (chrom, start, end, name) — lexicographic
    log.info(f'Sorting {len(rows)} total rows...')
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))

    log.info(f'Writing {args.output}...')
    with open(args.output, 'w') as out:
        for chrom, start, end, name, score, strand in rows:
            out.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n')

    log.info('Done.')


if __name__ == '__main__':
    main()
