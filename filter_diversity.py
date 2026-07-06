#!/usr/bin/env python3
"""
FILTER step: greedy max-min ("farthest point") pairwise-dissimilar selection
over a candidate manifest, by chain-CA RMSD. This is the same "diversity"
scoring target as the original pdz_pairwise_rmsd.py, generalized to consume
a manifest (any FIND step's output) and a configurable chain field/id
instead of a hard-coded directory glob and "chain A".

Writes rmsd_matrix.csv and a selection manifest into --out-dir, and (to keep
shell orchestration simple, matching the original script's convenience)
prints the selected structures' absolute pdb paths space-separated on
stdout.

Usage:
    python3 filter_diversity.py --in candidates.json --out-dir analysis \
        [--chain-field chain_domain] [--n-select 5]
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from manifest import read_manifest_stdin_or_path, write_manifest  # noqa: E402
from rmsd import (  # noqa: E402
    load_ca_objects, pairwise_rmsd_matrix, greedy_max_min_selection,
    require_single_chain_id, write_matrix_csv,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="in_path", required=True, help="Input candidate manifest ('-' for stdin)")
    parser.add_argument("--chain-field", default="chain_domain", help="Manifest field naming the chain id to align on (default: chain_domain)")
    parser.add_argument("--n-select", type=int, default=5, help="Number of dissimilar structures to pick")
    parser.add_argument("--out-dir", required=True, help="Directory to write rmsd_matrix.csv / selection.json")
    args = parser.parse_args()

    candidates = read_manifest_stdin_or_path(args.in_path)
    if not candidates:
        sys.exit("No candidates in input manifest")
    if len(candidates) < args.n_select:
        sys.exit(f"Only {len(candidates)} candidates but --n-select {args.n_select} requested")

    chain_id = require_single_chain_id(candidates, args.chain_field)

    os.makedirs(args.out_dir, exist_ok=True)

    ids = [c["id"] for c in candidates]
    if len(set(ids)) != len(ids):
        sys.exit("Duplicate candidate ids in manifest; run a best-score FILTER first to collapse duplicates")

    print(f"Loading chain-{chain_id} CA atoms for {len(candidates)} candidates ...", file=sys.stderr)
    named_paths = [(c["id"], c["pdb_path"]) for c in candidates]
    names = load_ca_objects(named_paths, chain_id)

    print(f"Computing pairwise cealign RMSD over {len(names)} structures ...", file=sys.stderr)

    def progress(done, total):
        print(f"  ...{done}/{total} pairs aligned", file=sys.stderr)

    matrix = pairwise_rmsd_matrix(names, progress=progress)

    matrix_csv = os.path.join(args.out_dir, "rmsd_matrix.csv")
    write_matrix_csv(matrix_csv, names, matrix)
    print(f"Wrote RMSD matrix to {matrix_csv}", file=sys.stderr)

    selected_idx = greedy_max_min_selection(matrix, args.n_select)
    selected = [candidates[i] for i in selected_idx]

    submatrix = {
        names[i]: {names[j]: matrix[i][j] for j in selected_idx}
        for i in selected_idx
    }
    selection_json = os.path.join(args.out_dir, "selection.json")
    write_manifest(selection_json, selected)
    with open(selection_json) as f:
        payload = json.load(f)
    payload["pairwise_rmsd"] = submatrix
    with open(selection_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Wrote selection to {selection_json}", file=sys.stderr)

    min_pairwise = min(
        matrix[i][j] for a, i in enumerate(selected_idx) for j in selected_idx[a + 1:]
    )
    print(f"Selected {len(selected)} structures (min pairwise RMSD = {min_pairwise:.2f} A):", file=sys.stderr)
    for c in selected:
        print(f"  {c['id']}  ({c['pdb_path']})", file=sys.stderr)

    print(" ".join(c["pdb_path"] for c in selected))


if __name__ == "__main__":
    main()
