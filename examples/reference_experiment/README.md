# reference_experiment

## Figure

A 2×1 comparison figure: a 1×5 "reference" montage (5 pairwise-dissimilar
structures from a plain directory) stacked above a 1×5 "experiment" montage
(5 pairwise-dissimilar structures among the best-scoring design per target
from a production campaign). Two independent SELECT → GENERATE → ASSEMBLE
branches combined into one final panel.

## Workflow

Two parallel branches, each FIND → FILTER → FILTER → GENERATE → ASSEMBLE,
joined by a final ASSEMBLE step:

**Reference branch** (`ref_*`):
1. **FIND** (`find_flat`) — enumerate a plain directory of PDBs.
2. **FILTER** (`filter_diversity`) — 5 pairwise-dissimilar structures by
   chain-domain Cα RMSD.
3. **GENERATE** (`generate_each` / `figures/pdz_figure.py`) — one render per
   selected structure.
4. **ASSEMBLE** (`assemble`, `no_scale: true`) — 1×5 montage, unscaled.

**Experiment branch** (`exp_*`):
1. **FIND** (`find_campaign`) — walk a production campaign layout, attaching
   confidence scores.
2. **FILTER** (`filter_best_score`) — keep the best-scoring design per
   target.
3. **FILTER** (`filter_diversity`) — 5 pairwise-dissimilar among those.
4. **GENERATE** (`generate_each` / `figures/pdz_figure.py`) — one render per
   selected structure.
5. **ASSEMBLE** (`assemble`, `no_scale: true`) — 1×5 montage, unscaled.

**Final**: `assemble` stacks the two montages 2×1 and is the *only* step
that scales to a fixed width — intermediate montages pass `--no-scale` so
repeated resizing never degrades image quality (`lib/montage.py`).

Available both as `reference_vs_experiment.yaml` (`run_pipeline.py`) and as
an equivalent hand-written `run_reference_vs_experiment.sh` wrapper taking
`--reference-dir`/`--prod-root`/`--out-dir`/etc. flags directly — a second
reference for how the same composition looks outside the YAML dispatcher.

## Run

```
cd examples/reference_experiment && python3 ../../run_pipeline.py reference_vs_experiment.yaml
```

or, the equivalent bash wrapper:

```
./run_reference_vs_experiment.sh [--reference-dir DIR] [--prod-root DIR] [--out-dir DIR]
```

No prerequisite steps.

## Files

| File | Purpose |
|---|---|
| `reference_vs_experiment.yaml` | The two-branch pipeline, run via `run_pipeline.py`. |
| `run_reference_vs_experiment.sh` | Equivalent hand-written bash wrapper, flag-driven. |
