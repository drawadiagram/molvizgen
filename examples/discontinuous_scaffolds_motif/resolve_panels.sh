#!/usr/bin/env bash
# Application-specific config: which 5 discontinuous-scaffolds MCSA-41 models
# go in the panel row (one representative per RESIDUE_ISLAND_COUNT in 2..6,
# picked from island_counts.csv among models that passed the pipeline's
# rmsd_threshold, per M0110_1c0p/M0157_1qh5/M0349_1e3v being requested
# explicitly and M0209_1lij/M0078_1al6 chosen as the best (lowest motif_rmsd)
# passing representatives for island counts 2 and 5), and where the
# completed production campaign data actually lives.
#
# Run from this directory; writes panel-spec JSONs (absolute paths) to
# $OUT_DIR/panel_specs, consumed by motif_panels_pipeline.yaml's render steps.
set -euo pipefail

PROD_ROOT="/home/mason/exdrive/rad/discontinuous_scaffolds/prod"
SCRIPTS_DIR="$PROD_ROOT/p1-p8r/discontinuous_scaffolds/scripts"
OUT_DIR="/home/mason/exdrive/rad/discontinuous_scaffolds/prod/molviz_motif_panels"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

python3 "$REPO_ROOT/examples/discontinuous_scaffolds_motif/find_best_fold.py" \
  --prod-root "$PROD_ROOT" \
  --design-json-dir "$SCRIPTS_DIR" \
  --reference-pdb-dir "$SCRIPTS_DIR/mcsa_41" \
  --island-counts-csv "$SCRIPTS_DIR/island_counts.csv" \
  --rmsd-threshold 1.5 \
  --models M0209_1lij M0110_1c0p M0349_1e3v M0078_1al6 M0157_1qh5 \
  --out-dir "$OUT_DIR/panel_specs"
