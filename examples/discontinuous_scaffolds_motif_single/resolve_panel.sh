#!/usr/bin/env bash
# Application-specific config: a single four-island-motif discontinuous-
# scaffolds model (M0097_1ctt — the best, i.e. lowest motif_rmsd, passing
# representative among the four-island models not already used in
# ../discontinuous_scaffolds_motif/'s five-panel row) for a close-up
# single-panel figure. Reuses the sibling example's find_best_fold.py
# lookup script rather than duplicating it.
#
# Run from this directory; writes the panel-spec JSON (absolute paths) to
# $OUT_DIR/panel_specs, consumed by motif_single_pipeline.yaml's render step.
set -euo pipefail

PROD_ROOT="/home/mason/exdrive/rad/discontinuous_scaffolds/prod"
SCRIPTS_DIR="$PROD_ROOT/p1-p8r/discontinuous_scaffolds/scripts"
OUT_DIR="/home/mason/exdrive/rad/discontinuous_scaffolds/prod/molviz_motif_single"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

python3 "$REPO_ROOT/examples/discontinuous_scaffolds_motif/find_best_fold.py" \
  --prod-root "$PROD_ROOT" \
  --design-json-dir "$SCRIPTS_DIR" \
  --reference-pdb-dir "$SCRIPTS_DIR/mcsa_41" \
  --island-counts-csv "$SCRIPTS_DIR/island_counts.csv" \
  --rmsd-threshold 1.5 \
  --models M0097_1ctt \
  --out-dir "$OUT_DIR/panel_specs"
