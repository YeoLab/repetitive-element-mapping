"""Convert a Gencode GTF to the 11-column parsed_ucsc_tableformat file."""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import open_maybe_gz, parse_gtf_attributes, setup_logger


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gtf', required=True, help='.gtf or .gtf.gz input file')
    p.add_argument('--output', required=True, help='destination file path')
    return p.parse_args()


def main():
    args = parse_args()
    log = setup_logger('generate_parsed_ucsc_tableformat')

    gtf_path = Path(args.gtf)
    if not gtf_path.exists():
        log.error(f"GTF file does not exist: {args.gtf}")
        sys.exit(1)

    # transcript data: keyed by (gene_id, transcript_id)
    # Each value: dict with chrom, strand, tx_start, tx_end, exons, cds_starts, cds_ends
    transcripts = {}
    records_read = 0

    log.info(f"Reading GTF: {args.gtf}")
    with open_maybe_gz(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attrs_field = parts[:9]
            if feature not in ('transcript', 'exon', 'CDS', 'stop_codon'):
                records_read += 1
                continue

            records_read += 1
            attrs = parse_gtf_attributes(attrs_field)
            gene_id = attrs.get('gene_id', '')
            transcript_id = attrs.get('transcript_id', '')
            key = (gene_id, transcript_id)

            start = int(start_s)
            end = int(end_s)
            # GTF is 1-based inclusive; convert to 0-based half-open
            start0 = start - 1
            end0 = end

            if key not in transcripts:
                transcripts[key] = {
                    'chrom': chrom,
                    'strand': strand,
                    'tx_start': None,
                    'tx_end': None,
                    'exons': [],
                    'cds_min': None,
                    'cds_max': None,
                }

            t = transcripts[key]

            if feature == 'transcript':
                t['tx_start'] = start0
                t['tx_end'] = end0
                t['chrom'] = chrom
                t['strand'] = strand
            elif feature == 'exon':
                t['exons'].append((start0, end0))
                if t['tx_start'] is None:
                    t['tx_start'] = start0
                    t['tx_end'] = end0
                else:
                    if start0 < t['tx_start']:
                        t['tx_start'] = start0
                    if end0 > t['tx_end']:
                        t['tx_end'] = end0
            elif feature in ('CDS', 'stop_codon'):
                if t['cds_min'] is None or start0 < t['cds_min']:
                    t['cds_min'] = start0
                if t['cds_max'] is None or end0 > t['cds_max']:
                    t['cds_max'] = end0

    log.info(f"Records read: {records_read}; transcripts found: {len(transcripts)}")

    header = '#ENSG\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds'

    log.info(f"Writing output: {args.output}")
    emitted = 0
    with open(args.output, 'w', encoding='utf-8') as out:
        out.write(header + '\n')
        for (gene_id, transcript_id) in sorted(transcripts.keys()):
            t = transcripts[(gene_id, transcript_id)]
            tx_start = t['tx_start']
            tx_end = t['tx_end']
            if tx_start is None:
                continue

            if t['cds_min'] is not None:
                cds_start = t['cds_min']
                cds_end = t['cds_max']
            else:
                # No CDS: UCSC convention is zero-length CDS at txStart
                cds_start = tx_start
                cds_end = tx_start

            exons = sorted(t['exons'])
            exon_count = len(exons)
            if exon_count == 0:
                exon_starts_str = f'{tx_start},'
                exon_ends_str = f'{tx_end},'
                exon_count = 1
            else:
                exon_starts_str = ','.join(str(s) for s, e in exons) + ','
                exon_ends_str = ','.join(str(e) for s, e in exons) + ','

            out.write(
                f'{gene_id}\t{transcript_id}\t{t["chrom"]}\t{t["strand"]}\t'
                f'{tx_start}\t{tx_end}\t{cds_start}\t{cds_end}\t'
                f'{exon_count}\t{exon_starts_str}\t{exon_ends_str}\n'
            )
            emitted += 1

    log.info(f"Transcripts emitted: {emitted}")


if __name__ == '__main__':
    main()
