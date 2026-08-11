#!/usr/bin/env python3
"""Parse RepBase 18.05 species-specific FASTA into repeat families.

The repeat portion of the bowtie2 index comes from
``RepBase18.05.fasta/species_specific/<species>_repbase*_v2.*``. Each record's
header is ``NAME_CLASS_TAXON``, and the index header is the NAME, uppercased.

Recovering NAME is the hard part, because both CLASS (``Simple_Repeat``) and
TAXON (``Homo_sapiens``) contain underscores. The rule here is:

1. **Strip TAXON** — the trailing binomial (``Genus_species``) if the last two
   tokens look like one, otherwise the single trailing token.
2. **Strip CLASS** — the *longest* matching suffix from
   ``refdata/repbase18.05-classes.txt``. Longest matters:
   ``CR1_HS_CR1_Homo_sapiens`` is NAME=``CR1_HS`` CLASS=``CR1``, so a
   left-to-right match returns the wrong name.
3. What remains is NAME.

Validated: this reproduces all 1,224 hg38 repeat-family names exactly, with
zero curated exceptions, and parses every non-bare mouse header.

**Mouse-specific hazards this module handles, which hg38 does not exhibit:**

- 11 mouse records have *bare* headers with no class or taxon (``>U1``, ``>UHG``).
  hg38 has none. A previous analysis concluded these families were absent from
  the mouse file because it searched for ``U1_``; they exist, just in a
  different format.
- Three of them (``U7``, ``U14``, ``U8``) duplicate a class-formatted record
  with a byte-identical sequence, so naive extraction emits the same family
  twice. `select_families` raises on any duplicate name rather than silently
  keeping one.
"""

import re
from pathlib import Path

REFDATA = Path(__file__).resolve().parents[3] / "refdata"
DROP_EXACT_FILE = REFDATA / "repbase18.05-drop-exact.txt"
CLASSES_FILE = REFDATA / "repbase18.05-classes.txt"

# Classes supplied by a dedicated source: SimpleRepeat kmers, the tRNA file and
# the RefSeq rRNA records. MamSINE1_tRNA_* is a tRNA-derived SINE, a real repeat
# family, and is kept.
DROP_CLASSES = ("Simple_Repeat", "tRNA", "rRNA")
KEEP_PREFIX = "MamSINE1_tRNA_"

_BINOMIAL = re.compile(r"^[A-Z][a-z]+_[a-z]+$")


def _load_list(path):
    return {
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def load_drop_exact(path=None):
    return _load_list(path or DROP_EXACT_FILE)


def load_classes(path=None):
    return _load_list(path or CLASSES_FILE)


def read_fasta(path):
    """Yield (header, sequence). Handles the mouse file's .fastq extension,
    which holds FASTA content despite the name."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def strip_taxon(header):
    """Remove the trailing binomial or single trailing token."""
    tokens = header.split("_")
    if len(tokens) >= 2 and _BINOMIAL.match("_".join(tokens[-2:])):
        return "_".join(tokens[:-2])
    return "_".join(tokens[:-1])


def family_name(header, classes):
    """Return the family NAME for a header, or None if it cannot be parsed.

    A bare header (no underscore) is its own name — the mouse file's U-snRNA
    records are stored that way.
    """
    if "_" not in header:
        return header
    body = strip_taxon(header)
    if not body:
        return None
    longest = None
    for cls in classes:
        if body.endswith("_" + cls) and (longest is None or len(cls) > len(longest)):
            longest = cls
    if longest is not None:
        return body[: -(len(longest) + 1)] or None
    return body


def keep_record(header, drop_exact):
    """Apply the family-selection rule (FIX-PLAN 4b)."""
    if header in drop_exact:
        return False
    if header.startswith(KEEP_PREFIX):
        return True
    return not any(f"_{cls}_" in header for cls in DROP_CLASSES)


def select_families(fasta_path, drop_exact=None, classes=None):
    """Return (families, unparsed).

    families: list of (name, header, sequence), name uppercased as the index
              stores it.
    unparsed: headers kept by the selection rule whose NAME could not be
              derived. Callers must treat a non-empty list as fatal rather than
              silently dropping families.

    Raises ValueError if two kept records resolve to the same family name.
    """
    drop_exact = load_drop_exact() if drop_exact is None else drop_exact
    classes = load_classes() if classes is None else classes

    families, unparsed, seen = [], [], {}
    for header, seq in read_fasta(fasta_path):
        if not keep_record(header, drop_exact):
            continue
        name = family_name(header, classes)
        if not name:
            unparsed.append(header)
            continue
        upper = name.upper()
        if upper in seen:
            raise ValueError(
                f"duplicate family {upper!r} from headers {seen[upper]!r} and "
                f"{header!r} — resolve in the drop-exact list before building"
            )
        seen[upper] = header
        families.append((upper, header, seq))
    return families, unparsed
