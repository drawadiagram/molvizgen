#!/usr/bin/env python3
"""
FILTER step: keep the top N candidates overall by --score-field, no
grouping. Complements filter_best_score.py (best-per-group): use this when
you want a fixed-size "top scorers" panel across an entire manifest instead
of one winner per group (e.g. picking 2 high-scoring designs for a figure,
regardless of which campaign they came from).

Usage:
    python3 filter_top_n.py --in candidates.json --out top.json \
        --score-field af2_plddt --n 2 [--mode max|min]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from manifest import read_manifest_stdin_or_path, write_manifest  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="in_path", required=True, help="Input candidate manifest ('-' for stdin)")
    parser.add_argument("--out", required=True, help="Path to write the filtered manifest JSON")
    parser.add_argument("--score-field", required=True, help="Field to rank candidates by (e.g. 'af2_plddt')")
    parser.add_argument("--n", type=int, required=True, help="Number of top candidates to keep")
    parser.add_argument("--mode", choices=["max", "min"], default="max", help="Keep the highest or lowest scorers (default: max)")
    args = parser.parse_args()

    candidates = read_manifest_stdin_or_path(args.in_path)
    scored = [c for c in candidates if c.get(args.score_field) is not None]
    dropped = len(candidates) - len(scored)
    if dropped:
        print(f"Dropped {dropped} candidates missing score field {args.score_field!r}", file=sys.stderr)

    scored.sort(key=lambda c: c[args.score_field], reverse=(args.mode == "max"))
    top = scored[:args.n]

    write_manifest(args.out, top)
    print(f"Kept top {len(top)} of {len(scored)} candidates by {args.mode} {args.score_field!r} -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
