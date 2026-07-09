# Examples

Each subdirectory here is a worked, runnable pipeline over one real (or
fixture-scale) data root, meant as both documentation of the YAML format
(see `run_pipeline.py`'s docstring and the main `CLAUDE.md`'s Architecture
section) and a template to copy for a new figure. Every example is run from
inside its own directory:

```
cd examples/<name> && python3 ../../run_pipeline.py <pipeline>.yaml
```

Some examples need an application-specific lookup step run first (they say
so below, and it's always a `./resolve_*.sh` or `./run_*.sh` wrapper in the
same directory) because the input they need — "which model", "which
production candidate" — isn't something a generic FIND/FILTER step can
answer; it depends on one particular external pipeline's directory layout
and CSV schema. Those lookup scripts write small JSON "spec" files that the
`run_pipeline.py` steps then read via a literal `input:`/`reference:`/
`design:` path — see each script's own docstring for the campaign-specific
knowledge it encodes.

All example YAMLs point at data roots under `/home/mason/exdrive/rad/...`
that are specific to this machine; treat the paths as illustrations of the
pattern, not portable defaults. Swap them (or use `--root`/`--out-root`, see
below) to point at your own data.

## heatmap

`rmsd_heatmap_pipeline.yaml` — the simplest example: FIND every `*.pdb` in a
directory, FILTER down to 5 pairwise-dissimilar structures by chain-domain
Cα RMSD (`filter_diversity`, backed by `lib/rmsd.py`), then plot the full
pairwise RMSD matrix as a heatmap PNG (`plot_heatmap`). Three steps, no
GENERATE/figure rendering at all — good first example to read.

## diverse_figures

`run_diverse_figures.sh` — a hand-written bash wrapper (not a
`run_pipeline.py` YAML) chaining `pdz_pairwise_rmsd.py` (the older,
directory-glob predecessor to `filter_diversity.py`) → `pdz_figure.py` (one
PyMOL render per selected structure) → `montage_figures.py`. Useful as a
second reference for how the scripts compose outside of the YAML dispatcher,
and for one-off runs: `./run_diverse_figures.sh [pdb-dir]`.

## reference_experiment

Two parallel SELECT → GENERATE → ASSEMBLE branches combined into one 2×1
figure: a "reference" branch (5 pairwise-dissimilar structures from a plain
directory) and an "experiment" branch (every candidate from a production
campaign, filtered to the best-scoring design per target via
`filter_best_score`, then to 5 pairwise-dissimilar among those). Each branch
renders its own 1×5 montage with `--no-scale`; only the final `assemble`
step scales to a fixed width, so repeated resizing never degrades image
quality (see `lib/montage.py`). Available both as
`reference_vs_experiment.yaml` (`run_pipeline.py`) and as an equivalent
hand-written `run_reference_vs_experiment.sh` wrapper that takes
`--reference-dir`/`--prod-root`/`--out-dir`/etc. flags directly.

## smallmol

`small_molecule_binder_comparison.yaml` — a five-panel figure over
`find_structures_smallmol.py`'s campaign layout: a left "problem" panel
(the target ligand alone, hot-spot atoms colored via
`ligand_hotspot_figure.py`) beside a 2×2 grid of designs (top row: 2
top-pLDDT designs from a "nonadaptive reference" campaign; bottom row: 2
top-pLDDT designs from the "prod" campaign), assembled with
`assemble_panel_layout` (a wide left panel + right grid, not expressible as
a uniform rows×cols montage). Demonstrates the `data:` block / `--root`
mechanism: the YAML declares `data.nonadaptive`/`data.prod` defaults, and the
same file can be re-pointed at a different campaign without editing it:

```
cd examples/smallmol && python3 ../../run_pipeline.py small_molecule_binder_comparison.yaml \
  --root prod=/path/to/prod/small_molecule_binding-1 \
  --root nonadaptive=/path/to/nonadaptive_reference/.../small_molecule_binding
```

## pdz_design_vs_template

Compares the pdzbinder campaign's single highest-confidence design against
the (usually much shorter-peptide) template structure its lineage was
originally generated from, both aligned on the conserved C-terminal "EPEA"
peptide motif (`lib/peptide_align.py`, Kabsch-fit on backbone atoms only —
no sequence/structure alignment needed since the correspondence is just
"last N residues"). Two views of the same pairing:
- `panel_4k6y` (`render_pair` / `aligned_pair_figure.py`): two separate,
  identically-oriented panels tiled side by side.
- `overlay_4k6y` (`render_overlay` / `aligned_overlay_figure.py`): one panel
  with both structures superposed (four colors, one per chain), plus solo
  exports of each structure alone since an overlay is easy to misread when
  the two cartoons occlude each other.

Prerequisite: `./resolve_template.sh` — resolves the winning candidate via
generic FIND (`find_structures_campaign.py`) + FILTER (`filter_top_n.py`),
then `resolve_template.py` (application-specific) looks up that candidate's
originating `_in` template PDB by stripping any `_subN` iterative-refinement
suffix from its `campaign_dir` (`find_structures_campaign.campaign_base`).
`run_design_vs_template.sh` runs both steps in sequence.

## discontinuous_scaffolds_motif

A five-panel row, one panel per `RESIDUE_ISLAND_COUNT` (2 through 6), each
showing a discontinuous-scaffolds design's protein (cartoon) Kabsch-aligned
onto a reference active site's ligand (licorice) with the design's motif
hot-spot atoms highlighted (spheres) — `motif_superposition_figure.py`,
built on `lib/rfd3_motif_select.py` + `lib/kabsch.py`.

Prerequisite: `./resolve_panels.sh` — `find_best_fold.py` (application
-specific) scans every `campaign_analysis.csv` under a production root for
each requested model's lowest-`motif_rmsd` row (mirroring the source
pipeline's own best-fold selection across however many redesign generations
a model went through) and resolves it to a Chai-1 CIF path, writing one
panel-spec JSON per model. `run_motif_panels.sh` runs both steps in
sequence.

## discontinuous_scaffolds_motif_single

A single close-up panel (`motif_single_pipeline.yaml`) for one
four-island-motif model (M0097_1ctt), reusing
`discontinuous_scaffolds_motif`'s `find_best_fold.py` rather than
duplicating it (see `resolve_panel.sh`). Exercises
`motif_superposition_figure.py`'s close-up options added specifically for
this figure without changing any prior default: `--ligand-representation
surface` (50%-transparent, so hot-spot spheres nested inside stay visible)
and `--orient-toward hotspot` (camera looks straight down on the best-fit
plane through just the motif atoms, rather than the whole-complex framing
the five-panel row uses). `run_motif_single.sh` runs `resolve_panel.sh` then
the pipeline.

## discontinuous_scaffolds_anchor_progression

A two-panel "before/after" row for a single input motif's redesign lineage:
the initial (root, failing) design next to the final (passing) redesign
generation, with the *anchor* residues — motif residues identified as
well-predicted enough to hold fixed during that redesign — highlighted in
both panels so the same physical residues are traceable across generations.
Ships in three renderer variants, each its own script/pipeline family (one
YAML per model, named `..._<model>.yaml`):
- `anchor_progression_pipeline*.yaml` (`motif_superposition_figure.py`):
  anchor residues as burnt-orange spheres/cartoon, rest of the hot spot
  purple; each panel re-orients itself independently.
- `anchor_progression_oriented_pipeline*.yaml`
  (`anchor_progression_oriented_figure.py`): both panels share one fixed
  camera orientation (looking down on the *initial* design's own motif
  plane) instead of re-orienting per panel; anchor residue green, rest of
  hot spot pink.
- `anchor_progression_zoned_pipeline*.yaml`
  (`anchor_progression_zoned_figure.py`): same shared orientation, plus the
  protein is colored by distance from the motif (near = opaque yellow
  cartoon, far = translucent white surface) instead of a flat color.

Prerequisite: `./resolve_progression.sh [MODEL]` (default `M0157_1qh5`;
`M0349_1e3v` and `M0739_1knp` are the other worked models with their own
YAMLs) — `find_anchor_progression.py` (application-specific) picks, among
every independent pipeline run that ever processed the requested model, the
one whose redesign lineage actually reaches a passing model (many don't
converge at all), then writes `<model>_initial.json`/`<model>_final.json`
panel specs. Every generation's spec references the *root* pipeline run's
design JSON, never a redesign's own — a redesign's contig/`select_fixed_atoms`
keys are renumbered relative to the true reference PDB, so only the root
contig's chai-sequence-position arithmetic (invariant across redesign
generations) can resolve a motif residue's position in a later generation's
fold. `run_anchor_progression*.sh` wrappers run the resolve step then a
specific pipeline variant.
