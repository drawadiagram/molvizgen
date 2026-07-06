#!/usr/bin/env python3
"""
PyMOL-backed chain-CA RMSD core, generalized from pdz_pairwise_rmsd.py so it
operates on an explicit list of (name, pdb_path) pairs and a caller-supplied
chain id, rather than a directory glob and a hard-coded "chain A".

This module has no directory or application-specific assumptions: it is pure
structural-comparison logic, reusable by any FILTER step that needs to reduce
a set of structures to a pairwise-dissimilar subset.
"""
import csv
import os
import sys

from pymol_scene import cmd


def load_ca_objects(named_paths, chain_id):
    """Load `chain_id`, CA atoms only, from each (name, pdb_path) pair.

    Returns the list of pymol object names actually loaded (same order as
    input), after verifying each has at least one atom post-filtering.
    """
    names = []
    for name, path in named_paths:
        cmd.load(path, name)
        cmd.remove(f"{name} and not chain {chain_id}")
        cmd.remove(f"{name} and not (polymer.protein and name CA)")
        if cmd.count_atoms(name) == 0:
            sys.exit(f"{name}: no chain-{chain_id} CA atoms found after filtering ({path})")
        names.append(name)
    return names


def pairwise_rmsd_matrix(names, progress=None):
    """Full symmetric RMSD matrix via cealign on the CA-only objects."""
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
            if progress and (done % 200 == 0 or done == total_pairs):
                progress(done, total_pairs)
    return matrix


def greedy_max_min_selection(matrix, n_select):
    """Greedy farthest-point selection: maximize the minimum pairwise RMSD
    among the chosen set, giving a pairwise-dissimilar subset. Returns a list
    of indices into `matrix`/`names`."""
    n = len(matrix)

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


def require_single_chain_id(candidates, chain_field):
    """Check that every candidate agrees on `chain_field`, and return that
    single agreed-upon chain id. Exits with an error message otherwise -
    every structure-comparison step here (filter_diversity.py,
    plot_rmsd_heatmap.py) needs one consistent chain id to align on."""
    chain_ids = {c.get(chain_field) for c in candidates}
    if len(chain_ids) != 1:
        sys.exit(f"Candidates disagree on {chain_field!r}: {chain_ids}; pass a manifest with a single consistent chain id")
    chain_id = chain_ids.pop()
    if chain_id is None:
        sys.exit(f"No {chain_field!r} field found on candidates")
    return chain_id


def write_matrix_csv(path, names, matrix):
    """Write a pairwise matrix to `path` as CSV: header row of `names`, then
    one row per name with its RMSD to every other name, formatted to 4
    decimal places."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + names)
        for name, row in zip(names, matrix):
            writer.writerow([name] + [f"{v:.4f}" for v in row])
