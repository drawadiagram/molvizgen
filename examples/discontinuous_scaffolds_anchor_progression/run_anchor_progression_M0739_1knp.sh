#!/usr/bin/env bash
# One-shot wrapper: resolve M0739_1knp's initial/final panel specs from the
# completed discontinuous-scaffolds campaign, then run the render+assemble
# pipeline. Equivalent to running `resolve_progression.sh M0739_1knp`
# followed by run_pipeline.py directly (see
# anchor_progression_pipeline_M0739_1knp.yaml's header comment) — same
# convention as run_anchor_progression.sh (the M0157_1qh5 counterpart) and
# run_anchor_progression_M0349_1e3v.sh (the M0349_1e3v counterpart).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

"$HERE/resolve_progression.sh" M0739_1knp
cd "$HERE" && python3 "$REPO_ROOT/run_pipeline.py" anchor_progression_pipeline_M0739_1knp.yaml
