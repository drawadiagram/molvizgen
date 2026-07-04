#!/usr/bin/env python3
"""
Compute pairwise structural RMSD between the PDZ domains (chain A, CA atoms)
of every PDB file in a directory, then pick a maximally-diverse subset via
greedy max-min ("farthest point") selection.

Usage:
    python3 pdz_pairwise_rmsd.py [--pdb-dir DIR] [--n-select 5] [--out-dir analysis]
"""
import argparse
import csv
import glob
import json
import os
import sys

import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd  # noqa: E402


def load_pdz_ca_objects(pdb_dir):
    """Load chain A (the PDZ domain), CA atoms only, from every *.pdb in pdb_dir."""
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdb_files:
        sys.exit(f"No .pdb files found in {pdb_dir}")

    names = []
    for path in pdb_files:
        name = os.path.splitext(os.path.basename(path))[0]
        cmd.load(path, name)
        cmd.remove(f"{name} and not chain A")
        cmd.remove(f"{name} and not (polymer.protein and name CA)")
        if cmd.count_atoms(name) == 0:
            sys.exit(f"{name}: no chain-A CA atoms found after filtering")
        names.append(name)
    return names


def pairwise_rmsd_matrix(names):
    """Full symmetric RMSD matrix via cealign on the CA-only chain-A objects."""
    n = len(names)
    matrix = [[0.0] * n for _ in range(n)]
    total_pairs = n * (n - 1) // 2
    done = 0
    for i in range(n):
        for j in range(i + 1, n):
            cmd.create("__mobile", names[i])
            result = cmd.cealign(names[j], "__mobile")
            cmd.delete("__mobile")
            rmsd = result["RMSD"]
            matrix[i][j] = matrix[j][i] = rmsd
            done += 1
            if done % 200 == 0 or done == total_pairs:
                print(f"  ...{done}/{total_pairs} pairs aligned", file=sys.stderr)
    return matrix


def greedy_max_min_selection(matrix, names, n_select):
    """Greedy farthest-point selection: maximize the minimum pairwise RMSD
    among the chosen set, giving a pairwise-dissimilar subset."""
    n = len(names)

    # Seed with the single most dissimilar pair.
    best_pair, best_rmsd = None, -1.0
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] > best_rmsd:
                best_rmsd = matrix[i][j]
                best_pair = (i, j)
    selected = list(best_pair)

    while len(selected) < n_select:
        best_candidate, best_min_dist = None, -1.0
        for k in range(n):
            if k in selected:
                continue
            min_dist_to_selected = min(matrix[k][s] for s in selected)
            if min_dist_to_selected > best_min_dist:
                best_min_dist = min_dist_to_selected
                best_candidate = k
        selected.append(best_candidate)

    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-dir", default=".", help="Directory containing the *.pdb files")
    parser.add_argument("--n-select", type=int, default=5, help="Number of dissimilar structures to pick")
    parser.add_argument("--out-dir", default="analysis", help="Directory to write matrix/selection output")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading PDZ domains (chain A, CA-only) from {args.pdb_dir} ...", file=sys.stderr)
    names = load_pdz_ca_objects(args.pdb_dir)
    print(f"Loaded {len(names)} structures. Computing pairwise cealign RMSD ...", file=sys.stderr)
    matrix = pairwise_rmsd_matrix(names)

    matrix_csv = os.path.join(args.out_dir, "rmsd_matrix.csv")
    with open(matrix_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + names)
        for name, row in zip(names, matrix):
            writer.writerow([name] + [f"{v:.4f}" for v in row])
    print(f"Wrote RMSD matrix to {matrix_csv}", file=sys.stderr)

    selected_idx = greedy_max_min_selection(matrix, names, args.n_select)
    selected_names = [names[i] for i in selected_idx]

    submatrix = {
        names[i]: {names[j]: matrix[i][j] for j in selected_idx}
        for i in selected_idx
    }
    selection_json = os.path.join(args.out_dir, "selected_five.json")
    with open(selection_json, "w") as f:
        json.dump(
            {"selected": selected_names, "pairwise_rmsd": submatrix},
            f,
            indent=2,
        )
    print(f"Wrote selection to {selection_json}", file=sys.stderr)

    min_pairwise = min(
        matrix[i][j] for a, i in enumerate(selected_idx) for j in selected_idx[a + 1:]
    )
    print(f"Selected {len(selected_names)} structures (min pairwise RMSD = {min_pairwise:.2f} A):", file=sys.stderr)
    for n in selected_names:
        print(f"  {n}", file=sys.stderr)

    # Space-separated on stdout for easy shell consumption.
    print(" ".join(selected_names))


if __name__ == "__main__":
    main()
