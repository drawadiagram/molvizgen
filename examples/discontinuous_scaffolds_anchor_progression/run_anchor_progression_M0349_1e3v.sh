#!/usr/bin/env bash
# One-shot wrapper: resolve M0349_1e3v's initial/final panel specs from the
# completed discontinuous-scaffolds campaign, then run the render+assemble
# pipeline. Equivalent to running `resolve_progression.sh M0349_1e3v`
# followed by run_pipeline.py directly (see
# anchor_progression_pipeline_M0349_1e3v.yaml's header comment) — same
# convention as run_anchor_progression.sh (the M0157_1qh5 counterpart).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

"$HERE/resolve_progression.sh" M0349_1e3v
cd "$HERE" && python3 "$REPO_ROOT/run_pipeline.py" anchor_progression_pipeline_M0349_1e3v.yaml
