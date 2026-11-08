#!/usr/bin/env python3
"""MetaX 2 adapter preparation: reject fractional DBs; make SE IDs unique.

MetaX groups mappings by read identifier. BAM-derived SE input may contain
two independent records with the same QNAME (e.g. mates whose pairing flags
were removed). Give each SE record a unique identifier without changing
sequence or record order. This is internal input, not a published profile.
PE input remains untouched: -i/-I carry the pairing.
"""

import argparse
import gzip
from pathlib import Path


def require_full_index(db):
    path = Path(db) / "metax_db.tsv.gz"
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#fractional\t"):
                if line.strip() != "#fractional\tfalse":
                    raise ValueError("MetaX adapter requires a full index; "
                                     "fractional databases are unsupported")
                return
    raise ValueError(f"{path}: missing #fractional metadata (expected MetaX 2)")


def unique_fasta_ids(source, destination):
    count = 0
    with open(source) as inp, open(destination, "w") as out:
        for line in inp:
            if line.startswith(">"):
                count += 1
                out.write(f">metax_se_{count}\n")
            else:
                if not count and line.strip():
                    raise ValueError(f"{source}: sequence before first header")
                out.write(line)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True)
    parser.add_argument("--single_end", choices=("true", "false"), required=True)
    parser.add_argument("--se", required=True)
    parser.add_argument("--out_se", required=True)
    args = parser.parse_args()
    try:
        require_full_index(args.db)
        if args.single_end == "true":
            unique_fasta_ids(args.se, args.out_se)
    except (OSError, ValueError) as error:
        parser.exit(1, f"[metax_prepare] ERROR: {error}\n")


if __name__ == "__main__":
    main()
