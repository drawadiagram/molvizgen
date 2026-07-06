#!/usr/bin/env bash
# One-shot wrapper: resolve the panel spec, then run the render pipeline.
# Same two-step composition as
# ../discontinuous_scaffolds_motif/run_motif_panels.sh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

"$HERE/resolve_panel.sh"
cd "$HERE" && python3 "$REPO_ROOT/run_pipeline.py" motif_single_pipeline.yaml
