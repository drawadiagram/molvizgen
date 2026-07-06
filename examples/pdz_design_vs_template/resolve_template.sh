#!/usr/bin/env bash
# Application-specific config: the pdzbinder production campaign this
# example draws from, and where to write the resolved analysis/template
# specs. FIND (find_structures_campaign.py) + FILTER (filter_top_n.py) do
# the general "every production complex, highest confidence_score first"
# work; resolve_template.py only adds the campaign-specific "_in" template
# lookup (see that script's docstring).
#
# Run from this directory; writes candidates_all.json/best.json under
# $OUT_DIR/analysis, and one <id>.json template spec (absolute paths) per
# selected candidate to $OUT_DIR/template_specs, consumed by
# design_vs_template_pipeline.yaml's render_pair step.
set -euo pipefail

PROD_ROOT="/home/mason/exdrive/rad/pdzbinder/prod"
OUT_DIR="/home/mason/exdrive/rad/pdzbinder/prod/molviz_design_vs_template"
N_BEST=1

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ANALYSIS_DIR="$OUT_DIR/analysis"
mkdir -p "$ANALYSIS_DIR"

echo "== FIND: every production complex across $PROD_ROOT ==" >&2
python3 "$REPO_ROOT/find_structures_campaign.py" \
    --prod-root "$PROD_ROOT" \
    --out "$ANALYSIS_DIR/candidates_all.json"

echo "== FILTER: top $N_BEST by confidence_score ==" >&2
python3 "$REPO_ROOT/filter_top_n.py" \
    --in "$ANALYSIS_DIR/candidates_all.json" \
    --score-field confidence_score --n "$N_BEST" --mode max \
    --out "$ANALYSIS_DIR/best.json"

echo "== Resolve each winner's design-template ('_in') structure ==" >&2
python3 "$(dirname "${BASH_SOURCE[0]}")/resolve_template.py" \
    --candidates "$ANALYSIS_DIR/best.json" \
    --prod-root "$PROD_ROOT" \
    --out-dir "$OUT_DIR/template_specs"
