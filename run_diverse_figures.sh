#!/usr/bin/env bash
# Pick 5 pairwise-dissimilar PDZ-domain structures (by chain-A CA RMSD),
# render a PyMOL figure for each, and cobble them into a single montage.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PDB_DIR="${1:-$(cd "$SCRIPT_DIR/../new_decoys" && pwd)}"
OUT_DIR="$(cd "$PDB_DIR" && pwd)"
ANALYSIS_DIR="$OUT_DIR/analysis"
FIGURES_DIR="$OUT_DIR/figures"

mkdir -p "$ANALYSIS_DIR" "$FIGURES_DIR"

echo "== Step 1/3: computing pairwise PDZ-domain RMSD and selecting 5 dissimilar structures ==" >&2
SELECTED=$(python3 "$SCRIPT_DIR/pdz_pairwise_rmsd.py" \
    --pdb-dir "$PDB_DIR" \
    --n-select 5 \
    --out-dir "$ANALYSIS_DIR")

echo "Selected: $SELECTED" >&2

echo "== Step 2/3: rendering PyMOL figures ==" >&2
FIGURE_PNGS=()
for id in $SELECTED; do
    pdb_path="$PDB_DIR/$id.pdb"
    out_png="$FIGURES_DIR/${id}_pdz_complex.png"
    echo "  rendering $id -> $out_png" >&2
    python3 "$SCRIPT_DIR/pdz_figure.py" "$pdb_path" "$out_png"
    FIGURE_PNGS+=("$out_png")
done

echo "== Step 3/3: building montage ==" >&2
MONTAGE_PNG="$FIGURES_DIR/montage.png"
python3 "$SCRIPT_DIR/montage_figures.py" "${FIGURE_PNGS[@]}" --out "$MONTAGE_PNG"

echo "" >&2
echo "Done. 5 figures written to $FIGURES_DIR:" >&2
for id in $SELECTED; do
    echo "  ${id}_pdz_complex.png"
done
echo "Montage: $MONTAGE_PNG" >&2
