#!/usr/bin/env bash
# One-shot wrapper: resolve the highest-confidence production design's
# template pairing from the completed pdzbinder campaign, then run the
# render+assemble pipeline. Equivalent to running resolve_template.sh
# followed by run_pipeline.py directly (see design_vs_template_pipeline.yaml's
# header comment) -- kept here as a second reference for how the two steps
# compose, same convention as
# examples/discontinuous_scaffolds_motif/run_motif_panels.sh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

"$HERE/resolve_template.sh"
cd "$HERE" && python3 "$REPO_ROOT/run_pipeline.py" design_vs_template_pipeline.yaml
