#!/usr/bin/env bash
# Application-specific config: which input motif and prod root to resolve
# a design-progression panel pair for. Takes the model name as an optional
# first argument (default M0157_1qh5, the RESIDUE_ISLAND_COUNT=6 example);
# M0349_1e3v (RESIDUE_ISLAND_COUNT=4) is the other model this example covers
# — see anchor_progression_pipeline_M0349_1e3v.yaml — chosen as the model
# with the most anchor residues identified at its root generation (2: A40,
# A100) among every model whose redesign lineage actually ends in a passing
# model (find_anchor_progression.py's docstring explains why a single model
# can have multiple independent pipeline runs, and why some redesign
# lineages never converge below rmsd_threshold at all).
#
# Run from this directory; writes panel-spec JSONs (absolute paths) to
# $OUT_DIR/panel_specs, consumed by anchor_progression_pipeline.yaml's (or
# anchor_progression_pipeline_M0349_1e3v.yaml's) render steps.
set -euo pipefail

PROD_ROOT="/home/mason/exdrive/rad/discontinuous_scaffolds/prod"
MODEL="${1:-M0157_1qh5}"
OUT_DIR="/home/mason/exdrive/rad/discontinuous_scaffolds/prod/molviz_anchor_progression"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

python3 "$REPO_ROOT/examples/discontinuous_scaffolds_anchor_progression/find_anchor_progression.py" \
  --prod-root "$PROD_ROOT" \
  --model "$MODEL" \
  --rmsd-threshold 1.5 \
  --out-dir "$OUT_DIR/panel_specs"
