#!/usr/bin/env bash
# One-shot wrapper: resolve the 5 models' panel specs from the completed
# discontinuous-scaffolds campaign, then run the render+assemble pipeline.
# Equivalent to running resolve_panels.sh followed by run_pipeline.py
# directly (see motif_panels_pipeline.yaml's header comment) — kept here as
# a second reference for how the two steps compose, same convention as
# examples/diverse_figures/run_diverse_figures.sh and
# examples/reference_experiment/run_reference_vs_experiment.sh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

"$HERE/resolve_panels.sh"
cd "$HERE" && python3 "$REPO_ROOT/run_pipeline.py" motif_panels_pipeline.yaml
