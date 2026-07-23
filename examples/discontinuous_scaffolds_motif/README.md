# discontinuous_scaffolds_motif

## Figure

A five-panel row, one panel per `RESIDUE_ISLAND_COUNT` (2 through 6), each
showing a discontinuous-scaffolds design's protein (cartoon, yellow)
Kabsch-aligned onto a reference active site's ligand (licorice, blue), with
the design's motif hot-spot atoms highlighted (spheres, purple).

## Workflow

One application-specific lookup step, then GENERATE (×5) → ASSEMBLE:

1. **Resolve** (`resolve_panels.sh`, prerequisite) — `find_best_fold.py`
   (application-specific) scans every `campaign_analysis.csv` under a
   production root for each requested model's lowest-`motif_rmsd` row
   (mirroring the source pipeline's own best-fold selection across however
   many redesign generations a model went through) and resolves it to a
   Chai-1 CIF path, writing one panel-spec JSON per model.
2. **GENERATE** (`render` / `figures/motif_superposition_figure.py`, ×5) —
   one panel per model: Kabsch-align the folded design's motif backbone
   atoms onto the reference active site (`lib/rfd3_motif_select.py` +
   `lib/kabsch.py`), then render protein cartoon + hot-spot spheres +
   reference ligand licorice.
3. **ASSEMBLE** (`assemble`) — `row`: tile the 5 panels into one 1×5 row.

## Run

```
./resolve_panels.sh
cd examples/discontinuous_scaffolds_motif && python3 ../../run_pipeline.py motif_panels_pipeline.yaml
```

or, both steps in sequence:

```
./run_motif_panels.sh
```

`motif_panels_pipeline_456.yaml` is a variant covering only islands 4-6.

## Files

| File | Purpose |
|---|---|
| `resolve_panels.sh` | Prerequisite: resolves 5 models' best folds to panel-spec JSONs. |
| `find_best_fold.py` | Application-specific: lowest-`motif_rmsd` fold lookup across a campaign. |
| `motif_panels_pipeline.yaml` | The five-panel (islands 2-6) row pipeline. |
| `motif_panels_pipeline_456.yaml` | Variant covering only islands 4-6. |
| `run_motif_panels.sh` | One-shot wrapper: `resolve_panels.sh` then `run_pipeline.py`. |
