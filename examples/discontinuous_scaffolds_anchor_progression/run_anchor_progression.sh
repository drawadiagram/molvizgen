#!/usr/bin/env bash
# One-shot wrapper: resolve M0157_1qh5's initial/final panel specs from the
# completed discontinuous-scaffolds campaign, then run the render+assemble
# pipeline. Equivalent to running resolve_progression.sh followed by
# run_pipeline.py directly (see anchor_progression_pipeline.yaml's header
# comment) — kept here as a second reference for how the two steps compose,
# same convention as examples/discontinuous_scaffolds_motif/run_motif_panels.sh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

"$HERE/resolve_progression.sh"
cd "$HERE" && python3 "$REPO_ROOT/run_pipeline.py" anchor_progression_pipeline.yaml
