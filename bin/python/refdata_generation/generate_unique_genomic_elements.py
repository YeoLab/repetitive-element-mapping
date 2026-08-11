"""Generate UniqueGenomicElements BED from RepeatMasker, tRNA, miRNA, and Gencode sources."""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from bin.python.refdata_generation._shared import open_maybe_gz, parse_gtf_attributes, setup_logger


# Families excluded from the unique-genome track: multi-copy rRNA and mitochondrial
# transcripts. Reads over these must be assigned to their repeat family, not to a
# "unique" genomic locus (cf. the RNA45S rRNA_extra_hash handling in
# parse_bowtie2_output_realtime_includemultifamily_{PE,SE}.pl). Derived from the hg38
# reference: these 5 families account for 582 of the 591 MASTER_FILELIST transcripts
# that the reference BED omits.
GENCODE_EXCLUDED_FAMILIES = frozenset({'RNA5S', 'RNA5-8S', 'MTTRNA', 'MTRNR1', 'MTRNR2'})

# Trailing gtRNAdb copy-number suffix, e.g. tRNA-Asn-GTT-2-3 -> tRNA-Asn-GTT
TRNA_COPY_SUFFIX = re.compile(r'-[0-9]+-[0-9]+$')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repeatmasker', required=True)
    p.add_argument('--trna')
    p.add_argument('--gff3')
    p.add_argument('--parsed-ucsc')
    p.add_argument('--master-filelist', required=True,
                   help='MASTER_FILELIST TSV; restricts the Gencode contribution to its '
                        'curated transcripts (col1=ENST, col4=family)')
    p.add_argument('--chrom-allowlist', required=True,
                   help='One scaffold name per line (or a .fai). Rows on any other '
                        'scaffold are dropped. Pin to the assembly patch release the '
                        'reference was built from (hg38 reference = GRCh38.p13).')
    p.add_argument('--assembly', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--flank', type=int, default=500)
    return p.parse_args()


def read_chrom_allowlist(path):
    """Read scaffold names from a plain list or a .fai (first whitespace field)."""
    names = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            names.add(line.split()[0])
    return names


def read_gencode_allowed_transcripts(path):
    """Transcript ids from MASTER_FILELIST whose family is not rRNA/mitochondrial."""
    allowed = set()
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4 or not parts[0].startswith('ENST'):
                continue
            if parts[3] in GENCODE_EXCLUDED_FAMILIES:
                continue
            allowed.add(parts[0])
    return allowed


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


def parse_trna_bed_rows(path):
    """Yield tRNA rows as (chrom, start0, end0, gene_id, family, strand).

    Column 5 carries the tRNA family (gene_id minus its copy-number suffix), not the
    GTF score -- the hg38 reference uses e.g. 'tRNA-Asn-GTT-2-3' / 'tRNA-Asn-GTT'.
    """
    for chrom, start0, end0, name, _score, strand in parse_gtf_bed_rows(path, name_field='gene_id'):
        yield (chrom, start0, end0, name, TRNA_COPY_SUFFIX.sub('', name), strand)


def ucsc_chrom(seqid):
    """Normalize a gff3 seqid to UCSC naming.

    miRBase mixes conventions. The GRCm39 mmu.gff3 (v23) is 1,164 rows of
    Ensembl-style seqids ('1', 'X') and 26 rows already written 'chr2', 'chr10'.
    Everything else in this repo -- the genome FASTAs, the chromosome
    allowlists, parsed_ucsc, the RepeatMasker tracks -- is UCSC-named, so an
    un-normalized seqid silently produces intervals that match no chromosome.
    """
    if seqid.startswith('chr'):
        return seqid
    return 'chrM' if seqid == 'MT' else 'chr' + seqid


def parse_gff3_mirna(path):
    """Parse miRNA_primary_transcript features from gff3. Yields (chrom, start0, end0, id, name, strand).

    Seqids are normalized through ucsc_chrom.
    """
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
            yield (ucsc_chrom(chrom), start0, end0, entry_id, entry_name, strand)


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

    allowed_chroms = read_chrom_allowlist(args.chrom_allowlist)
    log.info(f'Chromosome allowlist: {len(allowed_chroms)} scaffolds')
    allowed_transcripts = read_gencode_allowed_transcripts(args.master_filelist)
    log.info(f'Gencode allowlist: {len(allowed_transcripts)} curated transcripts '
             f'(families {sorted(GENCODE_EXCLUDED_FAMILIES)} excluded)')

    rows = []

    # RepeatMasker (required) — all instances, gene_id.
    # Simple repeats are deliberately NOT a source here: the hg38 reference BED contains
    # zero simple-repeat rows. They are represented as SimpleRepeat kmers in the bowtie2
    # index instead, so emitting them would double-count 1.05M loci as "unique genome".
    log.info('Reading RepeatMasker...')
    rm_rows = list(parse_gtf_bed_rows(args.repeatmasker, name_field='gene_id'))
    log.info(f'  {len(rm_rows)} RepeatMasker entries')
    rows.extend(rm_rows)

    # tRNA (optional) — gene_id in col4, family in col5
    if args.trna:
        log.info('Reading tRNA...')
        trna_rows = list(parse_trna_bed_rows(args.trna))
        log.info(f'  {len(trna_rows)} tRNA entries')
        rows.extend(trna_rows)
    else:
        log.warning(f'--trna not provided; tRNA entries omitted from {args.assembly} UniqueGenomicElements')

    # Gencode transcripts (optional) — transcript_id, score="-", actual strand;
    # restricted to the MASTER_FILELIST curated set
    if args.parsed_ucsc:
        log.info('Reading Gencode from parsed_ucsc...')
        gc_all = list(parse_parsed_ucsc(args.parsed_ucsc))
        gc_rows = [r for r in gc_all if r[3] in allowed_transcripts]
        log.info(f'  {len(gc_rows)} Gencode transcript entries kept of {len(gc_all)}')
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

    # Drop scaffolds absent from the pinned assembly release (newer _fix/_alt patches)
    before = len(rows)
    rows = [r for r in rows if r[0] in allowed_chroms]
    log.info(f'Chromosome allowlist dropped {before - len(rows)} of {before} rows')

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
