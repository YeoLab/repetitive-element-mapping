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
