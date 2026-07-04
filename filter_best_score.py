#!/usr/bin/env python3
"""
FILTER step: generic groupby-and-keep-best. Groups candidates in a manifest
by --key and keeps, from each group, the single record with the best
(max or min) --score-field. Candidates missing --score-field are dropped
with a warning (there is no score to compare).

This is the "the optimal designs are the best one from the campaign" filter:
e.g. collapsing several stage-variants of the same target down to the one
highest-confidence design, before a diversity filter picks a dissimilar
subset across targets.

Usage:
    python3 filter_best_score.py --in candidates.json --out best.json \
        --key id --score-field confidence_score [--mode max|min]
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
    parser.add_argument("--key", required=True, help="Field to group candidates by (e.g. 'id')")
    parser.add_argument("--score-field", required=True, help="Field to compare within each group (e.g. 'confidence_score')")
    parser.add_argument("--mode", choices=["max", "min"], default="max", help="Keep the max or min scoring candidate per group (default: max)")
    args = parser.parse_args()

    candidates = read_manifest_stdin_or_path(args.in_path)

    groups = {}
    dropped_no_key = 0
    dropped_no_score = 0
    for c in candidates:
        if args.key not in c:
            dropped_no_key += 1
            continue
        if c.get(args.score_field) is None:
            dropped_no_score += 1
            continue
        groups.setdefault(c[args.key], []).append(c)

    if dropped_no_key:
        print(f"Dropped {dropped_no_key} candidates missing key {args.key!r}", file=sys.stderr)
    if dropped_no_score:
        print(f"Dropped {dropped_no_score} candidates missing score field {args.score_field!r}", file=sys.stderr)

    pick = max if args.mode == "max" else min
    best = [pick(members, key=lambda c: c[args.score_field]) for members in groups.values()]
    best.sort(key=lambda c: c[args.key])

    write_manifest(args.out, best)
    print(f"Kept {len(best)} candidates (1 per {args.key}, {args.mode} {args.score_field}) -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
