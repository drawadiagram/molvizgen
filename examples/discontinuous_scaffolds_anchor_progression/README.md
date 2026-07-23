# discontinuous_scaffolds_anchor_progression

## Figure

A "before/after" row for a single input motif's redesign lineage: the
initial (root, failing) design next to the final (passing) redesign
generation, with the *anchor* residues — motif residues identified as
well-predicted enough to hold fixed during that redesign — highlighted in
both panels so the same physical residues are traceable across generations.

Ships in four renderer variants (each its own script), all over the same
resolved initial/final panel-spec pair, and one YAML per model per variant
(named `..._<model>.yaml`):

- **`anchor_progression_pipeline*.yaml`**
  (`figures/motif_superposition_figure.py`) — anchor residues as
  burnt-orange spheres/cartoon, rest of the hot spot purple; each panel
  re-orients itself independently.
- **`anchor_progression_oriented_pipeline*.yaml`**
  (`figures/anchor_progression_oriented_figure.py`) — both panels share one
  fixed camera orientation (looking down on the *initial* design's own
  motif plane) instead of re-orienting per panel; anchor residue green,
  rest of hot spot pink.
- **`anchor_progression_zoned_pipeline*.yaml`**
  (`figures/anchor_progression_zoned_figure.py`) — same shared orientation,
  plus the protein is colored by distance from the motif (near = opaque
  yellow cartoon, far = translucent white surface) instead of a flat color.
- **`anchor_progression_motif_pipeline*.yaml`**
  (`figures/anchor_progression_oriented_figure.py --motif-only`) — same
  shared orientation and hot-spot/anchor coloring as the oriented variant,
  but drops the whole-protein cartoon entirely: only the motif residues (as
  sticks) plus the reference ligand are shown, camera zoomed tight to just
  that content. Unlike the other three variants, there is no `assemble`
  step — the two panels are kept as standalone images rather than montaged
  into a row.

## Workflow

One application-specific lookup step per model, then GENERATE (×2, one
`render_pair` step) → ASSEMBLE (skipped for the motif-only variant):

1. **Resolve** (`resolve_progression.sh [MODEL]`, prerequisite, default
   `M0157_1qh5`) — `find_anchor_progression.py` (application-specific)
   picks, among every independent pipeline run that ever processed the
   requested model, the one whose redesign lineage actually reaches a
   passing model (many don't converge at all), then writes
   `<model>_initial.json`/`<model>_final.json` panel specs. Every
   generation's spec references the *root* pipeline run's design JSON,
   never a redesign's own — a redesign's contig/`select_fixed_atoms` keys
   are renumbered relative to the true reference PDB, so only the root
   contig's chai-sequence-position arithmetic (invariant across redesign
   generations) can resolve a motif residue's position in a later
   generation's fold.
2. **GENERATE** (`render` ×2, or `render_pair` for the oriented/zoned/
   motif-only variants) — `panel_initial`/`panel_final` (or `panel_pair`):
   render the root and passing-redesign generations with the same
   hot-spot/anchor highlighting.
3. **ASSEMBLE** (`assemble`) — `row`: tile the two panels side by side.
   Skipped by the motif-only variant, which has no wider complex to frame.

`M0157_1qh5` (RESIDUE_ISLAND_COUNT=6) is the default model; `M0349_1e3v`
(RESIDUE_ISLAND_COUNT=4) and `M0739_1knp` are the other worked models with
dedicated `run_anchor_progression_<model>.sh` wrappers for the base
(`anchor_progression_pipeline*`) variant — `M0110_1c0p`, `M0711_2esd`, and
`M0738_1o98` have YAMLs for the oriented/motif-only variants but no
dedicated wrapper script (run `resolve_progression.sh <model>` then
`run_pipeline.py <variant>_pipeline_<model>.yaml` directly).

## Run

```
./resolve_progression.sh [MODEL]          # default M0157_1qh5
cd examples/discontinuous_scaffolds_anchor_progression && python3 ../../run_pipeline.py anchor_progression_pipeline.yaml
# or: anchor_progression_oriented_pipeline_<model>.yaml
# or: anchor_progression_zoned_pipeline_<model>.yaml
# or: anchor_progression_motif_pipeline_<model>.yaml
```

or, for the base variant's three worked models:

```
./run_anchor_progression.sh              # M0157_1qh5
./run_anchor_progression_M0349_1e3v.sh
./run_anchor_progression_M0739_1knp.sh
```

## Files

| File | Purpose |
|---|---|
| `resolve_progression.sh` | Prerequisite: resolves a model's initial/final panel specs. Takes model name as an optional first argument. |
| `find_anchor_progression.py` | Application-specific: picks the passing redesign lineage for a model. |
| `anchor_progression_pipeline*.yaml` | Base variant: independent per-panel orientation, spheres/cartoon coloring. |
| `anchor_progression_oriented_pipeline*.yaml` | Shared fixed camera orientation across both panels. |
| `anchor_progression_zoned_pipeline*.yaml` | Shared orientation + distance-from-motif protein coloring. |
| `anchor_progression_motif_pipeline*.yaml` | Shared orientation, motif-only close-up (no whole-protein cartoon, no assemble step). |
| `run_anchor_progression*.sh` | One-shot wrappers: `resolve_progression.sh` then a specific base-variant pipeline. |
