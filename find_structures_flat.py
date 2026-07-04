#!/usr/bin/env python3
"""
FIND step: enumerate every *.pdb file in a flat directory and emit a
candidate manifest (see lib/manifest.py). Pure discovery — no scoring.

Usage:
    python3 find_structures_flat.py --dir DIR --out candidates.json \
        [--chain-domain A] [--chain-peptide B] [--glob '*.pdb']
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from manifest import write_manifest  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", required=True, help="Directory to search for structure files")
    parser.add_argument("--glob", default="*.pdb", help="Glob pattern within --dir (default: *.pdb)")
    parser.add_argument("--chain-domain", default="A", help="Chain id to record for the domain (default: A)")
    parser.add_argument("--chain-peptide", default="B", help="Chain id to record for the peptide (default: B)")
    parser.add_argument("--out", required=True, help="Path to write the candidate manifest JSON")
    args = parser.parse_args()

    pdb_dir = os.path.abspath(args.dir)
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, args.glob)))
    if not pdb_files:
        sys.exit(f"No files matching {args.glob!r} found in {pdb_dir}")

    candidates = []
    for path in pdb_files:
        cid = os.path.splitext(os.path.basename(path))[0]
        candidates.append({
            "id": cid,
            "pdb_path": path,
            "chain_domain": args.chain_domain,
            "chain_peptide": args.chain_peptide,
        })

    write_manifest(args.out, candidates)
    print(f"Found {len(candidates)} candidates in {pdb_dir} -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
