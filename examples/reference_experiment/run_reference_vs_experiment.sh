#!/usr/bin/env bash
# Reference vs. experiment workflow: build a 1x5 montage of the 5
# pairwise-dissimilar structures found in a reference directory (this
# application's existing SELECT|GENERATE|ASSEMBLE pipeline), build a second
# 1x5 montage of the 5 pairwise-dissimilar, best-scoring structures found in
# a production campaign, then ASSEMBLE the two into a single 2x1 figure.
#
# Only the final ASSEMBLE step scales to a fixed width — the two
# intermediate 1x5 panels are assembled at native resolution (--no-scale).
#
# Usage:
#   run_reference_vs_experiment.sh \
#       --reference-dir DIR --prod-root DIR --out-dir DIR \
#       [--n-select 5] [--groups 'p1-p16,p17-p32'] [--target-width-in 6.0]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

REFERENCE_DIR="$(pwd)"
PROD_ROOT="/home/mason/exdrive/rad/pdzbinder/prod"
OUT_DIR="$(pwd)"
N_SELECT=5
PROD_GROUPS=""
TARGET_WIDTH_IN=6.0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference-dir) REFERENCE_DIR="$2"; shift 2 ;;
        --prod-root) PROD_ROOT="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --n-select) N_SELECT="$2"; shift 2 ;;
        --groups) PROD_GROUPS="$2"; shift 2 ;;
        --target-width-in) TARGET_WIDTH_IN="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

REFERENCE_DIR="$(cd "$REFERENCE_DIR" && pwd)"
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

ANALYSIS_DIR="$OUT_DIR/analysis"
FIGURES_DIR="$OUT_DIR/figures"
mkdir -p "$ANALYSIS_DIR/reference" "$ANALYSIS_DIR/experiment" "$FIGURES_DIR/reference" "$FIGURES_DIR/experiment"

# --- Branch 1: reference (this directory's existing pipeline) ------------
echo "== Reference branch: SELECT (FIND + FILTER) ==" >&2
python3 "$REPO_ROOT/find_structures_flat.py" \
    --dir "$REFERENCE_DIR" \
    --out "$ANALYSIS_DIR/reference/candidates.json"

REF_SELECTED=$(python3 "$REPO_ROOT/filter_diversity.py" \
    --in "$ANALYSIS_DIR/reference/candidates.json" \
    --chain-field chain_domain \
    --n-select "$N_SELECT" \
    --out-dir "$ANALYSIS_DIR/reference")

echo "== Reference branch: GENERATE ==" >&2
REF_PNGS=()
for pdb_path in $REF_SELECTED; do
    id="$(basename "$pdb_path" .pdb)"
    out_png="$FIGURES_DIR/reference/${id}_pdz_complex.png"
    echo "  rendering $id -> $out_png" >&2
    python3 "$REPO_ROOT/figures/pdz_figure.py" "$pdb_path" "$out_png"
    REF_PNGS+=("$out_png")
done

echo "== Reference branch: ASSEMBLE (unscaled) ==" >&2
REF_MONTAGE="$FIGURES_DIR/reference_montage.png"
python3 "$REPO_ROOT/montage_figures.py" "${REF_PNGS[@]}" \
    --rows 1 --cols "$N_SELECT" --out "$REF_MONTAGE" --no-scale

# --- Branch 2: experiment (production campaign) ---------------------------
echo "== Experiment branch: SELECT (FIND + FILTER best-score + FILTER diversity) ==" >&2
FIND_ARGS=(--prod-root "$PROD_ROOT" --out "$ANALYSIS_DIR/experiment/candidates_all.json")
if [[ -n "$PROD_GROUPS" ]]; then
    FIND_ARGS+=(--groups "$PROD_GROUPS")
fi
python3 "$REPO_ROOT/find_structures_campaign.py" "${FIND_ARGS[@]}"

python3 "$REPO_ROOT/filter_best_score.py" \
    --in "$ANALYSIS_DIR/experiment/candidates_all.json" \
    --out "$ANALYSIS_DIR/experiment/candidates_best.json" \
    --key id --score-field confidence_score --mode max

EXP_SELECTED=$(python3 "$REPO_ROOT/filter_diversity.py" \
    --in "$ANALYSIS_DIR/experiment/candidates_best.json" \
    --chain-field chain_domain \
    --n-select "$N_SELECT" \
    --out-dir "$ANALYSIS_DIR/experiment")

# find_structures_campaign.py already normalized these to standard chain
# A (domain) / B (peptide), so the existing, unmodified pdz_figure.py (which
# defaults to exactly that convention) renders them with no extra flags —
# the chain-naming difference was resolved upstream at FIND time, not here.
echo "== Experiment branch: GENERATE ==" >&2
EXP_PNGS=()
for pdb_path in $EXP_SELECTED; do
    id="$(basename "$pdb_path" .pdb)"
    out_png="$FIGURES_DIR/experiment/${id}_pdz_complex.png"
    echo "  rendering $id -> $out_png" >&2
    python3 "$REPO_ROOT/figures/pdz_figure.py" "$pdb_path" "$out_png"
    EXP_PNGS+=("$out_png")
done

echo "== Experiment branch: ASSEMBLE (unscaled) ==" >&2
EXP_MONTAGE="$FIGURES_DIR/experiment_montage.png"
python3 "$REPO_ROOT/montage_figures.py" "${EXP_PNGS[@]}" \
    --rows 1 --cols "$N_SELECT" --out "$EXP_MONTAGE" --no-scale

# --- Final: ASSEMBLE the two panels, scaling only here ---------------------
echo "== Final ASSEMBLE: reference (top) vs. experiment (bottom), scaled to ${TARGET_WIDTH_IN}in ==" >&2
FINAL_PNG="$FIGURES_DIR/reference_vs_experiment.png"
python3 "$REPO_ROOT/montage_figures.py" "$REF_MONTAGE" "$EXP_MONTAGE" \
    --rows 2 --cols 1 --out "$FINAL_PNG" --target-width-in "$TARGET_WIDTH_IN"

echo "" >&2
echo "Done. Final figure: $FINAL_PNG" >&2
