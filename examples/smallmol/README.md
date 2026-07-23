# smallmol

## Figure

A five-panel small-molecule-binder comparison: a left "problem" panel (the
target ligand alone, hot-spot atoms red vs. the rest yellow, per an
RFDiffusion3 binder-design spec) beside a 2×2 grid of designs — top row the
best- and median-AF2-pLDDT designs from a "nonadaptive reference" campaign,
bottom row the same best+median pick from the "prod" campaign (one
cherry-picked winner plus one representative-of-the-distribution design per
campaign, not two winners). Each design panel carries the problem panel's
red-hot-spot/yellow-rest ligand coloring through onto its bound peptide, so
the same atoms read as targeted vs. not in every panel.

## Workflow

One SELECT/GENERATE branch for the problem panel, two parallel FIND → FILTER
→ GENERATE branches for the design grid, joined by one ASSEMBLE step:

1. **GENERATE** (`render` / `figures/ligand_hotspot_figure.py`) —
   `problem_panel`: the target ligand alone, split into hot-spot/rest atoms.
2. **FIND** (`find_smallmol`, once per campaign) — discover completed
   designs from an ImpressBasePipeline campaign layout, pairing each AF2
   fold-prediction with its preceding FastRelax protein+ligand complex.
3. **FILTER** (`filter_top_n --select top_and_median`, once per campaign) —
   a fixed best-scorer-plus-median-scorer pick by AF2 pLDDT (`--n`
   ignored).
4. **GENERATE** (`generate_each` / `figures/generate_figure.py`, once per
   campaign) — one render per picked design, colored by an RFDiffusion3
   binder-design spec's hot-spot split (`--peptide-design-json`) so the
   bound peptide carries the problem panel's coloring.
5. **ASSEMBLE** (`assemble_panel_layout`) — the problem panel (1/3 width) next
   to the 2×2 design grid (2/3 width), all four design panels tiled in one
   `build_montage()` call so cell size/spacing is consistent across both
   rows.

Demonstrates the `data:` block / `--root` mechanism: the YAML declares
`data.nonadaptive`/`data.prod` defaults, and the same file can be re-pointed
at a different campaign without editing it.

## Run

```
cd examples/smallmol && python3 ../../run_pipeline.py small_molecule_binder_comparison.yaml
```

Points at campaign `-2` by default. Re-point at a different campaign
without editing the YAML:

```
cd examples/smallmol && python3 ../../run_pipeline.py small_molecule_binder_comparison.yaml \
  --root prod=/path/to/prod/small_molecule_binding-1 \
  --root nonadaptive=/path/to/nonadaptive_reference/.../small_molecule_binding
```

No prerequisite steps.

## Files

| File | Purpose |
|---|---|
| `small_molecule_binder_comparison.yaml` | The whole five-panel pipeline. |
