"""Generate a combined FASTA and bowtie2 index from Gencode, RepeatMasker, tRNA, miRNA, and custom sequences."""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pybedtools

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import (
    assert_writable, faidx_if_missing, open_maybe_gz, parse_gtf_attributes, setup_logger,
)


# Non-coding RNA transcript types to include from Gencode (verified against hg38 reference)
DEFAULT_TRANSCRIPT_TYPES = {
    'snRNA', 'misc_RNA', 'snoRNA', 'rRNA_pseudogene', 'scaRNA',
    'rRNA', 'Mt_tRNA', 'Mt_rRNA', 'vaultRNA',
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gtf', required=True)
    p.add_argument('--parsed-ucsc', required=True)
    p.add_argument('--repeatmasker', required=True)
    p.add_argument('--simplerepeats', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--trna')
    p.add_argument('--gff3')
    p.add_argument('--custom-fasta', action='append', default=[])
    p.add_argument('--output-dir', required=True)
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--transcript-types', nargs='*', default=None,
                   help='Gencode transcript types to include. Default: non-coding RNA types.')
    return p.parse_args()


def _is_simple_repeat_gene_id(gid):
    return bool(re.match(r'^\(.+\)n$', gid))


def _simple_repeat_to_fasta_name(gid):
    """Convert '(AT)n' → 'AT_SimpleRepeat'."""
    inner = gid[1:-2]  # strip outer ( and )n
    return inner.upper() + '_SimpleRepeat'


def read_gtf_transcript_types(gtf_path):
    """Return dict: transcript_id → transcript_type from GTF transcript features."""
    tid_to_type = {}
    with open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'transcript':
                continue
            attrs = parse_gtf_attributes(parts[8])
            tid = attrs.get('transcript_id', '')
            tt = attrs.get('transcript_type', attrs.get('transcript_biotype', ''))
            if tid and tt:
                tid_to_type[tid] = tt
    return tid_to_type


def read_parsed_ucsc(path):
    """Return dict: transcript_id → (chrom, strand, txStart, txEnd, [(exon_start, exon_end), ...])."""
    transcripts = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            _, tid, chrom, strand, txstart, txend, _, _, _, exon_starts_s, exon_ends_s = parts[:11]
            starts = [int(x) for x in exon_starts_s.rstrip(',').split(',') if x]
            ends = [int(x) for x in exon_ends_s.rstrip(',').split(',') if x]
            exons = sorted(zip(starts, ends))
            transcripts[tid] = (chrom, strand, int(txstart), int(txend), exons)
    return transcripts


def build_bed12(transcripts):
    """Build BED12 string for all transcripts."""
    lines = []
    for tid, (chrom, strand, txstart, txend, exons) in transcripts.items():
        if not exons:
            continue
        block_count = len(exons)
        block_sizes = ','.join(str(e - s) for s, e in exons) + ','
        block_starts = ','.join(str(s - txstart) for s, e in exons) + ','
        lines.append(
            f'{chrom}\t{txstart}\t{txend}\t{tid}\t0\t{strand}\t{txstart}\t{txend}\t0\t{block_count}\t{block_sizes}\t{block_starts}'
        )
    return '\n'.join(lines) + '\n'


def extract_sequences_bed12(bed12_str, genome_fa, log, valid_chroms=None):
    """Use pybedtools to extract spliced sequences from BED12. Returns fasta string."""
    if valid_chroms is not None:
        lines = [l for l in bed12_str.splitlines(keepends=True) if not l or l.split('\t')[0] in valid_chroms]
        bed12_str = ''.join(lines)
    if not bed12_str.strip():
        return ''
    bt = pybedtools.BedTool(bed12_str, from_string=True)
    result = bt.sequence(fi=str(genome_fa), s=True, name=True, split=True)
    return open(result.seqfn).read()


def load_fai_chroms(fasta_path):
    """Return set of chromosome names present in the .fai index."""
    fai = str(fasta_path) + '.fai'
    chroms = set()
    with open(fai) as fh:
        for line in fh:
            chroms.add(line.split('\t')[0])
    return chroms


def extract_sequences_bed6(bed_rows, genome_fa, log, valid_chroms=None):
    """Extract sequences from 6-col BED rows. Returns fasta string."""
    if valid_chroms is not None:
        filtered = [r for r in bed_rows if r[0] in valid_chroms]
        skipped = len(bed_rows) - len(filtered)
        if skipped:
            log.debug(f'Skipped {skipped} BED rows on chromosomes not in FASTA')
        bed_rows = filtered
    bed_str = ''.join(
        f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n'
        for chrom, start, end, name, score, strand in bed_rows
    )
    if not bed_str.strip():
        return ''
    bt = pybedtools.BedTool(bed_str, from_string=True)
    result = bt.sequence(fi=str(genome_fa), s=True, name=True)
    return open(result.seqfn).read()


def read_repeatmasker(path, log):
    """Parse repeatmasker TSV. Returns (non_simple_rows, simple_rows).

    non_simple_rows: list of (chrom, start0, end0, UPPERCASE_gene_id, score, strand)
                     deduplicated by UPPERCASE_gene_id (first occurrence)
    simple_rows: list of (chrom, start0, end0, PATTERN_SimpleRepeat_name, score, strand)
                 deduplicated by (pattern)n gene_id (first occurrence)
    """
    seen_nonsimple = {}  # uppercase name → row
    seen_simple = {}     # (pattern)n original → row

    with open_maybe_gz(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attrs = parts[:9]
            if feature not in ('exon',):
                continue
            attr = parse_gtf_attributes(attrs)
            gid = attr.get('gene_id', '')
            if not gid:
                continue
            start0 = int(start_s) - 1
            end0 = int(end_s)
            score_val = score if score != '.' else '0'

            if _is_simple_repeat_gene_id(gid):
                if gid not in seen_simple:
                    name = _simple_repeat_to_fasta_name(gid)
                    seen_simple[gid] = (chrom, start0, end0, name, score_val, strand)
            else:
                up = gid.upper()
                if up not in seen_nonsimple:
                    seen_nonsimple[up] = (chrom, start0, end0, up, score_val, strand)

    return list(seen_nonsimple.values()), list(seen_simple.values())


def read_trna(path, log):
    """Parse tRNA TSV. Returns list of (chrom, start0, end0, gene_id, score, strand)."""
    rows = {}
    with open_maybe_gz(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attrs = parts[:9]
            if feature != 'exon':
                continue
            attr = parse_gtf_attributes(attrs)
            gid = attr.get('gene_id', '')
            if not gid or gid in rows:
                continue
            start0 = int(start_s) - 1
            end0 = int(end_s)
            score_val = score if score != '.' else '0'
            rows[gid] = (chrom, start0, end0, gid, score_val, strand)
    return list(rows.values())


def read_mirna_gff3(path, log):
    """Parse miRNA gff3. Returns list of (chrom, start0, end0, name, score, strand)."""
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attrs = parts[:9]
            # Only use primary miRNA entries (miRNA_primary_transcript) or miRNA
            if feature not in ('miRNA', 'miRNA_primary_transcript'):
                continue
            # Parse Name= from gff3 attributes
            m = re.search(r'Name=([^;]+)', attrs)
            if not m:
                continue
            name = m.group(1)
            if name in rows:
                continue
            start0 = int(start_s) - 1
            end0 = int(end_s)
            score_val = score if score != '.' else '0'
            rows[name] = (chrom, start0, end0, name, score_val, strand)
    return list(rows.values())


def read_custom_fasta(paths):
    """Concatenate custom FASTA files verbatim."""
    out = []
    for p in paths:
        with open(p) as fh:
            out.append(fh.read())
    return ''.join(out)


def count_fasta_seqs(fasta_str):
    return fasta_str.count('>')


def main():
    args = parse_args()
    log = setup_logger('generate_bowtie2_index')

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    assert_writable(output_dir, allowed_prefixes=[
        'examples/inputs/mm10', 'examples/inputs/mm39',
        'examples/inputs/hg38/bowtie2_index', '/tmp',
    ])

    genome_fa = Path(args.fasta)
    if not genome_fa.exists():
        log.error(f'Genome FASTA not found: {genome_fa}')
        sys.exit(1)
    faidx_if_missing(genome_fa)
    valid_chroms = load_fai_chroms(genome_fa)
    log.info(f'Genome has {len(valid_chroms)} chromosomes/scaffolds in .fai')

    out_fa = output_dir / (args.output_prefix + '.fa')

    # ── 1. Gencode transcripts ──────────────────────────────────────────
    transcript_types_filter = set(args.transcript_types) if args.transcript_types is not None else DEFAULT_TRANSCRIPT_TYPES
    log.info(f'Reading GTF for transcript types (filter: {len(transcript_types_filter)} types)...')
    tid_to_type = read_gtf_transcript_types(args.gtf)
    log.info(f'  {len(tid_to_type)} transcripts found in GTF')

    log.info('Reading parsed_ucsc_tableformat...')
    all_transcripts = read_parsed_ucsc(args.parsed_ucsc)
    transcripts = {tid: v for tid, v in all_transcripts.items()
                   if tid_to_type.get(tid, '') in transcript_types_filter}
    log.info(f'  {len(transcripts)} / {len(all_transcripts)} transcripts after type filter')
    bed12_str = build_bed12(transcripts)
    log.info('Extracting Gencode transcript sequences (BED12+split)...')
    gencode_fa = extract_sequences_bed12(bed12_str, genome_fa, log, valid_chroms)
    n_gencode = count_fasta_seqs(gencode_fa)
    log.info(f'  Extracted {n_gencode} / {len(transcripts)} transcript sequences')

    # ── 2 & 3. RepeatMasker (non-simple + simple) ───────────────────────
    log.info('Reading repeatmasker...')
    rm_rows, simple_rows = read_repeatmasker(args.repeatmasker, log)
    log.info(f'  Non-simple repeat elements: {len(rm_rows)}')
    log.info(f'  Simple repeat elements: {len(simple_rows)}')

    log.info('Extracting RepeatMasker sequences...')
    rm_fa = extract_sequences_bed6(rm_rows, genome_fa, log, valid_chroms)
    simple_fa = extract_sequences_bed6(simple_rows, genome_fa, log, valid_chroms)
    n_rm = count_fasta_seqs(rm_fa)
    n_simple = count_fasta_seqs(simple_fa)
    log.info(f'  Extracted {n_rm} repeat, {n_simple} simple-repeat sequences')

    # ── 4. tRNA (optional) ──────────────────────────────────────────────
    trna_fa = ''
    if args.trna:
        log.info('Reading tRNA...')
        trna_rows = read_trna(args.trna, log)
        log.info(f'  {len(trna_rows)} tRNA entries')
        trna_fa = extract_sequences_bed6(trna_rows, genome_fa, log, valid_chroms)
        log.info(f'  Extracted {count_fasta_seqs(trna_fa)} tRNA sequences')
    else:
        log.warning('--trna not provided; tRNA entries omitted')

    # ── 5. miRNA (optional) ─────────────────────────────────────────────
    mirna_fa = ''
    if args.gff3:
        log.info('Reading miRNA gff3...')
        mirna_rows = read_mirna_gff3(args.gff3, log)
        log.info(f'  {len(mirna_rows)} miRNA entries')
        mirna_fa = extract_sequences_bed6(mirna_rows, genome_fa, log, valid_chroms)
        log.info(f'  Extracted {count_fasta_seqs(mirna_fa)} miRNA sequences')
    else:
        log.warning('--gff3 not provided; miRNA entries omitted')

    # ── 6. Custom FASTA ─────────────────────────────────────────────────
    custom_fa = ''
    if args.custom_fasta:
        log.info(f'Reading {len(args.custom_fasta)} custom FASTA file(s)...')
        custom_fa = read_custom_fasta(args.custom_fasta)
        log.info(f'  {count_fasta_seqs(custom_fa)} custom sequences')

    # ── Write combined FASTA ────────────────────────────────────────────
    log.info(f'Writing combined FASTA: {out_fa}')
    with open(out_fa, 'w') as out:
        out.write(gencode_fa)
        out.write(rm_fa)
        out.write(simple_fa)
        out.write(trna_fa)
        out.write(mirna_fa)
        out.write(custom_fa)

    total_seqs = count_fasta_seqs(
        gencode_fa + rm_fa + simple_fa + trna_fa + mirna_fa + custom_fa
    )
    log.info(f'Total sequences written: {total_seqs}')

    # ── Missing ID report ───────────────────────────────────────────────
    expected = len(transcripts) + len(rm_rows) + len(simple_rows)
    if args.trna:
        expected += len(trna_rows)
    if args.gff3:
        expected += len(mirna_rows)
    extracted = n_gencode + n_rm + n_simple + count_fasta_seqs(trna_fa) + count_fasta_seqs(mirna_fa)
    missing = expected - extracted
    if expected > 0 and missing / expected > 0.01:
        report_path = output_dir / 'missing_ids_report.txt'
        log.warning(f'{missing}/{expected} IDs missing (>{1}%). Writing report: {report_path}')
        with open(report_path, 'w') as rpt:
            rpt.write(f'Missing {missing} of {expected} expected sequences\n')
    else:
        log.info(f'Missing sequences: {missing} / {expected}')

    # ── Build bowtie2 index ──────────────────────────────────────────────
    log.info('Building bowtie2 index...')
    index_prefix = str(output_dir / args.output_prefix)
    subprocess.run(['bowtie2-build', str(out_fa), index_prefix], check=True)

    log.info('Verifying bowtie2 index...')
    subprocess.run(
        ['bowtie2-inspect', '--summary', index_prefix],
        check=True, stdout=subprocess.DEVNULL,
    )
    log.info('Done.')


if __name__ == '__main__':
    main()
