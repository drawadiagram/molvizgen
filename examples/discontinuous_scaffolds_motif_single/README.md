# discontinuous_scaffolds_motif_single

## Figure

A single close-up panel for one four-island-motif discontinuous-scaffolds
design (M0097_1ctt): protein cartoon (yellow) and motif hot-spot atoms
(spheres, purple) as in `../discontinuous_scaffolds_motif`, but the
reference ligand rendered as a 50%-transparent surface (blue) instead of
licorice, and the camera looking straight down onto the best-fit plane of
the motif hot-spot atoms instead of the whole-complex framing the five-panel
row uses.

## Workflow

One application-specific lookup step, then a single GENERATE step (no
ASSEMBLE — there's only one panel):

1. **Resolve** (`resolve_panel.sh`, prerequisite) — reuses
   `../discontinuous_scaffolds_motif/find_best_fold.py` rather than
   duplicating it, restricted to one model, writing one panel-spec JSON.
2. **GENERATE** (`render` / `figures/motif_superposition_figure.py`) —
   `panel`: same Kabsch-alignment core as the five-panel row, exercising
   `--ligand-representation surface` (so the nested hot-spot spheres stay
   visible) and `--orient-toward hotspot` (camera faces the motif's own
   best-fit plane rather than the whole complex) — both added without
   changing either flag's default, so the five-panel pipeline keeps
   rendering identically. The camera still pans/zooms to fit the whole
   reference+design complex first; the close-up framing comes from a safe
   2D content-crop afterward (`--trim-pad`), never a tighter 3D camera move.

## Run

```
./resolve_panel.sh
cd examples/discontinuous_scaffolds_motif_single && python3 ../../run_pipeline.py motif_single_pipeline.yaml
```

or, both steps in sequence:

```
./run_motif_single.sh
```

## Files

| File | Purpose |
|---|---|
| `resolve_panel.sh` | Prerequisite: resolves M0097_1ctt's best fold to a panel-spec JSON. |
| `motif_single_pipeline.yaml` | The single close-up panel pipeline. |
| `run_motif_single.sh` | One-shot wrapper: `resolve_panel.sh` then `run_pipeline.py`. |
