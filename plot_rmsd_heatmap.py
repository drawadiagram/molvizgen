#!/usr/bin/env python3
"""
Render the pairwise chain-CA RMSD values computed by lib/rmsd.py as a
heatmap PNG.

Two input modes:

  --matrix rmsd_matrix.csv
      Plot a matrix already written by filter_diversity.py or
      pdz_pairwise_rmsd.py (first row/column are structure names, cells are
      RMSD in Angstrom). No PyMOL needed.

  --in candidates.json [--chain-field chain_domain] [--out-matrix FILE]
      Compute the matrix from scratch via lib/rmsd.py (load_ca_objects +
      pairwise_rmsd_matrix), like filter_diversity.py does, then plot it.
      Pass --out-matrix to also save the computed matrix as CSV.

Usage:
    python3 plot_rmsd_heatmap.py --matrix analysis/rmsd_matrix.csv --out heatmap.png
    python3 plot_rmsd_heatmap.py --in candidates.json --out heatmap.png --out-matrix matrix.csv
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from manifest import read_manifest_stdin_or_path, write_manifest  # noqa: E402,F401

# Sequential "blue" ramp, light -> dark (references/palette.md, dataviz skill):
# lightest step reads as "near zero", darkest as "most dissimilar".
BLUE_SEQUENTIAL = [
    "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
    "#256abf", "#1c5cab", "#104281", "#0d366b",
]
RMSD_CMAP = LinearSegmentedColormap.from_list("rmsd_blue", BLUE_SEQUENTIAL)

MUTED_INK = "#898781"
PRIMARY_INK = "#0b0b0b"


def read_matrix_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    names = rows[0][1:]
    matrix = [[float(v) for v in row[1:]] for row in rows[1:]]
    return names, matrix


def compute_matrix(in_path, chain_field):
    from rmsd import load_ca_objects, pairwise_rmsd_matrix  # noqa: E402

    candidates = read_manifest_stdin_or_path(in_path)
    if not candidates:
        sys.exit("No candidates in input manifest")

    chain_ids = {c.get(chain_field) for c in candidates}
    if len(chain_ids) != 1:
        sys.exit(f"Candidates disagree on {chain_field!r}: {chain_ids}; pass a manifest with a single consistent chain id")
    chain_id = chain_ids.pop()
    if chain_id is None:
        sys.exit(f"No {chain_field!r} field found on candidates")

    named_paths = [(c["id"], c["pdb_path"]) for c in candidates]
    print(f"Loading chain-{chain_id} CA atoms for {len(candidates)} candidates ...", file=sys.stderr)
    names = load_ca_objects(named_paths, chain_id)

    print(f"Computing pairwise cealign RMSD over {len(names)} structures ...", file=sys.stderr)

    def progress(done, total):
        print(f"  ...{done}/{total} pairs aligned", file=sys.stderr)

    matrix = pairwise_rmsd_matrix(names, progress=progress)
    return names, matrix


def write_matrix_csv(path, names, matrix):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + names)
        for name, row in zip(names, matrix):
            writer.writerow([name] + [f"{v:.4f}" for v in row])


def plot_heatmap(names, matrix, out_path, title, dpi, annotate):
    arr = np.array(matrix)
    n = len(names)

    fig_size = max(4.0, min(0.35 * n + 1.5, 20.0))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")

    im = ax.imshow(arr, cmap=RMSD_CMAP, vmin=0.0)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=max(5, 10 - n // 15), color=MUTED_INK)
    ax.set_yticklabels(names, fontsize=max(5, 10 - n // 15), color=MUTED_INK)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if annotate is None:
        annotate = n <= 20
    if annotate:
        vmax = arr.max() if arr.size else 1.0
        for i in range(n):
            for j in range(n):
                value = arr[i, j]
                text_color = "#fcfcfb" if vmax and value > 0.6 * vmax else PRIMARY_INK
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=max(5, 9 - n // 10), color=text_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Cα RMSD (Å)", color=PRIMARY_INK)
    cbar.ax.tick_params(colors=MUTED_INK)
    cbar.outline.set_visible(False)

    ax.set_title(title, color=PRIMARY_INK, fontsize=13, pad=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--matrix", help="Path to a precomputed rmsd_matrix.csv")
    source.add_argument("--in", dest="in_path", help="Candidate manifest to compute the matrix from ('-' for stdin)")
    parser.add_argument("--chain-field", default="chain_domain", help="Manifest field naming the chain id to align on, used with --in (default: chain_domain)")
    parser.add_argument("--out-matrix", help="With --in, also write the computed matrix to this CSV path")
    parser.add_argument("--out", required=True, help="Output heatmap PNG path")
    parser.add_argument("--title", default="Pairwise Cα RMSD", help="Plot title")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--annotate", action=argparse.BooleanOptionalAction, default=None, help="Print RMSD values in each cell (default: on for <=20 structures)")
    args = parser.parse_args()

    if args.matrix:
        names, matrix = read_matrix_csv(args.matrix)
    else:
        names, matrix = compute_matrix(args.in_path, args.chain_field)
        if args.out_matrix:
            write_matrix_csv(args.out_matrix, names, matrix)
            print(f"Wrote RMSD matrix to {args.out_matrix}", file=sys.stderr)

    plot_heatmap(names, matrix, args.out, args.title, args.dpi, args.annotate)
    print(f"Wrote {args.out} ({len(names)}x{len(names)} structures)")


if __name__ == "__main__":
    main()
