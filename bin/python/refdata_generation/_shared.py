"""Shared utility functions for refdata_generation scripts."""

import gzip
import logging
import os
import re
import subprocess
import sys
from pathlib import Path


def open_maybe_gz(path):
    """Open a .gtf or .gtf.gz file transparently, returning a text stream."""
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def parse_gtf_attributes(attrs_field):
    """Parse a GTF column-9 attribute string into a dict, handling 'key "value";' format."""
    result = {}
    for match in re.finditer(r'(\w+)\s+"([^"]*)"', attrs_field):
        result[match.group(1)] = match.group(2)
    return result


def setup_logger(name):
    """Return a logger that writes ISO-timestamped messages to stderr."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s [%(name)s]: %(message)s',
                                               datefmt='%Y-%m-%dT%H:%M:%S'))
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    return logger


TRNA_HEADER_RE = re.compile(r'^(\S+).*?(\S+):(\d+)-(\d+)\s+\(([+-])\)')


def read_gtrnadb_fasta(path):
    """Parse a gtRNAdb genomic tRNA FASTA.

    Returns list of (name, seq, chrom, start0, end0, strand), where name is the
    header's first token with the leading '<Genus>_<species>_' stripped, e.g.
    'Homo_sapiens_tRNA-Ala-AGC-1-1' → 'tRNA-Ala-AGC-1-1'. Coordinates come from
    the trailing 'chr6:28795964-28796035 (-)' field of the same header.
    """
    entries = []
    header, chunks = None, []

    def flush():
        if header is None:
            return
        m = TRNA_HEADER_RE.match(header)
        if not m:
            raise ValueError(f'Unparseable gtRNAdb header: {header!r}')
        raw, chrom, a, b, strand = m.groups()
        idx = raw.find('tRNA-')
        if idx < 0:
            raise ValueError(f'gtRNAdb name has no tRNA- component: {raw!r}')
        lo, hi = sorted((int(a), int(b)))
        entries.append((raw[idx:], ''.join(chunks).upper(), chrom, lo - 1, hi, strand))

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                flush()
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line)
    flush()
    return entries


def assert_writable(out_path, allowed_prefixes=None):
    """Raise PermissionError if out_path is not under any of allowed_prefixes."""
    if allowed_prefixes is None:
        allowed_prefixes = ['examples/inputs/mm10', 'examples/inputs/mm39', '/tmp']
    out_path = str(Path(out_path).resolve())
    for prefix in allowed_prefixes:
        abs_prefix = str(Path(prefix).resolve()) if not prefix.startswith('/') else prefix
        if out_path.startswith(abs_prefix):
            return
    raise PermissionError(
        f"Output path '{out_path}' is not within allowed prefixes: {allowed_prefixes}"
    )


def faidx_if_missing(fasta):
    """If <fasta>.fai is missing, invoke samtools faidx as subprocess."""
    fai = str(fasta) + '.fai'
    if not os.path.exists(fai):
        subprocess.run(['samtools', 'faidx', str(fasta)], check=True)
