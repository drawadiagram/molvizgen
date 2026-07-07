#!/usr/bin/env bash
# Application-specific config: which input motif and prod root to resolve
# a design-progression panel pair for. M0157_1qh5 has multiple independent
# pipeline runs (see find_anchor_progression.py's docstring); this auto-
# selects the one whose redesign lineage ends in a passing model
# (motif_rmsd < rmsd_threshold), which turns out to be the p11-p31r batch's
# disco_p14_0 (initial) -> disco_p14_R (final, passing) lineage.
#
# Run from this directory; writes panel-spec JSONs (absolute paths) to
# $OUT_DIR/panel_specs, consumed by anchor_progression_pipeline.yaml's
# render steps.
set -euo pipefail

PROD_ROOT="/home/mason/exdrive/rad/discontinuous_scaffolds/prod"
MODEL="M0157_1qh5"
OUT_DIR="/home/mason/exdrive/rad/discontinuous_scaffolds/prod/molviz_anchor_progression"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

python3 "$REPO_ROOT/examples/discontinuous_scaffolds_anchor_progression/find_anchor_progression.py" \
  --prod-root "$PROD_ROOT" \
  --model "$MODEL" \
  --rmsd-threshold 1.5 \
  --out-dir "$OUT_DIR/panel_specs"
