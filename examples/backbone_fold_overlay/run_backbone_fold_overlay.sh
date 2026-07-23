#!/usr/bin/env bash
# One-shot wrapper for backbone_fold_overlay_figure.py: render a de novo
# RFDiffusion3 backbone against the downstream AlphaFold2/ColabFold fold it
# led to, in a production small_molecule_binding-2 campaign (p2, task73 rfd3
# -> task78 alphafold -- chosen by surveying every rfd3->alphafold lineage
# across p1-p8 and picking the best structural agreement: RMSD 0.41 A over
# 96/99 residues against the fold's best-matching repeat window, since the
# fold prediction is of the whole C5-symmetric assembly the backbone seeds,
# folded as one long single chain -- see lib/oligomer_align.py).
#
# Produces all seven images this example's figure script supports: the
# backbone alone, the fold alone (cartoon), the two aligned overlays
# (backbone-ribbon/fold-cartoon and the swapped backbone-cartoon/fold-alt),
# the fold alone in --alt-representation (surface), the fold alone in
# ribbon, and the fold alone in licorice colored by secondary structure.
#
# Usage: ./run_backbone_fold_overlay.sh [backbone.cif.gz] [fold.pdb] [out_dir]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

CAMPAIGN_ROOT="/home/mason/exdrive/rad/smallmol/prod/small_molecule_binding-2/p2"
BACKBONE_PDB="${1:-$CAMPAIGN_ROOT/73_rfd3/out/ALR_binder_design_partial_0_model_0.cif.gz}"
FOLD_PDB="${2:-$CAMPAIGN_ROOT/78_alphafold/out/binder_unrelaxed_rank_001_alphafold2_model_3_seed_999.pdb}"
OUT_DIR="${3:-/home/mason/exdrive/rad/smallmol/backbone_fold_overlay/final}"

mkdir -p "$OUT_DIR"

python3 "$REPO_ROOT/backbone_fold_overlay_figure.py" \
    "$BACKBONE_PDB" "$FOLD_PDB" \
    "$OUT_DIR/backbone_alone.png" \
    "$OUT_DIR/fold_alone.png" \
    "$OUT_DIR/overlap.png" \
    "$OUT_DIR/overlap_swapped.png" \
    --out-fold-alt "$OUT_DIR/fold_surface_alone.png" \
    --out-fold-ribbon "$OUT_DIR/fold_ribbon_alone.png" \
    --out-fold-ss "$OUT_DIR/fold_licorice_by_ss.png"

echo "" >&2
echo "Done. Images written to $OUT_DIR:" >&2
for f in backbone_alone fold_alone overlap overlap_swapped fold_surface_alone fold_ribbon_alone fold_licorice_by_ss; do
    echo "  $f.png" >&2
done
