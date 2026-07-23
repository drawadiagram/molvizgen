# molvizgen

Command-line tools for taking a set of protein structure (PDB) files through
a **FIND → FILTER → GENERATE → ASSEMBLE** pipeline (with an optional **PLOT**
step for RMSD heatmaps): enumerate candidate structures, filter them down by
score or structural diversity, render publication-style PyMOL figures, and
assemble the results into montages — built from small, composable scripts,
driven either directly or by a YAML pipeline description.

## Usage

```bash
# FIND: enumerate a directory of PDB files into a manifest
python3 find_structures_flat.py --dir /path/to/pdbs --chain-domain A --chain-peptide B --out candidates.json

# FILTER: keep a diverse, maximally-dissimilar subset (by chain-CA RMSD)
python3 filter_diversity.py --in candidates.json --chain-field chain_domain --n-select 5 --out-dir analysis

# FILTER: keep only the best-scoring candidate per group
python3 filter_best_score.py --in candidates.json --key target --score-field confidence_score --mode max --out best.json

# FILTER: keep the top N candidates overall, no grouping
python3 filter_top_n.py --in candidates.json --score-field af2_plddt --n 10 --mode max --out top10.json

# FILTER: or keep a fixed best-scorer + median-scorer pair instead of top N
python3 filter_top_n.py --in candidates.json --score-field af2_plddt --select top_and_median --out top2.json

# PLOT: render a saved RMSD matrix as a heatmap
python3 plot_rmsd_heatmap.py --matrix analysis/rmsd_matrix.csv --title "Pairwise Cα RMSD" --out heatmap.png

# GENERATE: render one auto-oriented figure for a single structure
python3 pdz_figure.py my_structure.pdb figures/my_structure.png

# ASSEMBLE: tile rendered figures into a montage
python3 montage_figures.py figures/*.png --rows 2 --cols 3 --out montage.png

# Or run an entire FIND -> FILTER -> GENERATE -> ASSEMBLE pipeline from one YAML file
cd examples/heatmap && python3 ../../run_pipeline.py rmsd_heatmap_pipeline.yaml

# Or explore the pipeline interactively in a notebook
marimo edit rmsd_heatmap_notebook.py
```

## Concepts

**Manifest.** The one data structure every step agrees on: a JSON object
`{"candidates": [{"id": ..., "pdb_path": ..., "chain_domain": ..., ...}, ...]}`.
`pdb_path` is always absolute. A FIND step produces a manifest from a
directory scan; a FILTER step reads a manifest and writes a smaller or
reordered one; a GENERATE step reads a manifest and produces images. Because
every candidate carries its own absolute path and metadata, no step ever has
to reconstruct a file location from an id and a directory convention.

**Consistent framing.** Figures are rendered with an automatic camera: PCA
the domain chain's Cα atoms to find its long axis, rotate that to vertical,
then twist around that axis so the peptide/ligand chain faces the camera.
Every structure ends up framed the same way regardless of its original
coordinate frame, so a montage of different structures reads as one
consistent figure set rather than a grab-bag of arbitrary orientations.
`ligand_hotspot_figure.py` applies the same long-axis-to-vertical idea to a
lone ligand, since there's no second chain to twist toward the camera. This
orientation math (`lib/geometry.py`'s PCA/rotation helpers, `lib/orient.py`'s
scene-transform functions) and the shared PyMOL rendering boilerplate
(`lib/pymol_scene.py`: session bootstrap, material settings, ray-trace-and-save)
live in `lib/`, so every GENERATE script is a thin, domain-specific layer over
the same shared core rather than reimplementing it.

**Content-trimmed montages.** Every GENERATE script ray-traces onto a
fixed-size canvas with its own zoom buffer, so raw PNGs carry inconsistent
background margin — most visible when panels have different aspect ratios
(e.g. a tall thin ligand next to a roughly square protein cartoon).
`lib/imgtrim.py` crops each source image to its content bounding box (plus a
small pad) before tiling; `lib/montage.py`'s `build_montage` does this by
default (`--no-trim` to skip), shared by both `montage_figures.py` and
`assemble_panel_layout.py` rather than each carrying its own tiling logic.

## The pipeline

A pipeline is a YAML file: a required `figure_name` label, an ordered list of
named steps (each with a `kind` and `args`), and — separately — the data it
operates on. Every `kind` belongs to one of five stages: **FIND** (enumerate
a directory or campaign layout into a manifest), **FILTER** (reduce a
manifest by score or structural diversity), **PLOT** (render a saved RMSD
matrix as a heatmap), **GENERATE** (render PyMOL figures, per-candidate or
one-off), and **ASSEMBLE** (tile rendered figures into a montage or panel
layout) — see [Step kinds](#step-kinds) below for the full list. The example
below chains FIND → FILTER → PLOT:

```yaml
figure_name: rmsd_heatmap_demo
out_root: analysis

steps:
  - name: find
    kind: find_flat
    args:
      dir: ../new_decoys          # a sibling directory of PDB files

  - name: select
    kind: filter_diversity
    args:
      in: "${find.manifest}"      # reference an earlier step's output
      chain_field: chain_domain
      n_select: 5

  - name: heatmap
    kind: plot_heatmap
    args:
      matrix: "${select.matrix}"
      title: "Pairwise Cα RMSD"
      out: rmsd_heatmap.png
```

Run it with:
```
cd examples/heatmap && python3 ../../run_pipeline.py rmsd_heatmap_pipeline.yaml
```

`${step_name.field}` substitutes an earlier step's declared output — it must
be the *entire* value of a YAML key (not embedded in a larger string), and
resolves to whatever type that field holds (a path, or a list of paths for
things like a set of rendered PNGs). Each step gets its own output
subdirectory, `<out_root>/<figure_name>/<step name>/`.

### Data roots and output location

A pipeline YAML is meant to describe the **transformation** only — which
steps, filters, colors, layout — not hardcode the absolute paths of the data
it happens to have been developed against. Two mechanisms keep data out of
the transformation:

**Input data** — an optional top-level `data:` block declares named
directory defaults, referenced from any step's `args` as `${data.NAME}`:

```yaml
data:
  prod: /path/to/prod/experiment-1

steps:
  - name: prod_find
    kind: find_smallmol
    args:
      campaign_root: "${data.prod}"                     # whole-value
      # or, joined with a literal subpath in the same string:
      # input: "${data.prod}/p1_in/design.json"
```

Unlike `${step_name.field}`, `${data.NAME}` may be embedded inside a larger
string (e.g. joined with a subpath), since a data root is always a plain
string. Override or supply a root at invocation time with `--root
[NAME=]PATH` (repeatable; bare `PATH` with no `NAME=` sets the root named
`default`) — this is what makes one YAML reusable against a different
experiment directory with no edits:

```
python3 run_pipeline.py pipeline.yaml --root prod=/path/to/prod/experiment-2
```

**Output location** — `figure_name` (required) is a descriptive label used
as the output subfolder name; `out_root` (optional, top-level, defaults to
cwd `.`) is the base directory it's created under. Override `out_root` per
invocation with `--out-root PATH`, independent of any `data:` roots:

```
python3 run_pipeline.py pipeline.yaml --root prod=/path/to/experiment-2 --out-root /path/to/scratch
```

An undeclared, unsupplied `${data.NAME}` is a hard error at the step that
references it, with a message naming the missing root and how to supply it.

### Step kinds

| `kind` | Purpose | Key `args` | Declares |
|---|---|---|---|
| `find_flat` | FIND: enumerate every PDB file in a directory | `dir`, `glob` (default `*.pdb`), `chain_domain`, `chain_peptide` | `manifest` |
| `find_campaign` | FIND: recursively discover structures from a production campaign layout, attach confidence scores, normalize non-standard chain ids | `prod_root`, `groups`, `raw_chain_domain`, `raw_chain_peptide`, `normalized_dir` | `manifest` |
| `find_smallmol` | FIND: discover completed small-molecule-binder designs from an ImpressBasePipeline campaign layout, pairing each AF2 fold-prediction with its preceding FastRelax protein+ligand complex | `campaign_root`, `campaign_label`, `chain_domain`, `chain_peptide` | `manifest` |
| `filter_best_score` | FILTER: keep the best-scoring candidate per group (e.g. one design per target) | `in`, `key`, `score_field`, `mode` (`max`\|`min`) | `manifest` |
| `filter_top_n` | FILTER: keep the top N candidates overall by a score field, no grouping (or, `select: top_and_median`, a fixed best+median pick) | `in`, `score_field`, `n`, `mode` (`max`\|`min`), `select` (`top_n`\|`top_and_median`) | `manifest` |
| `filter_diversity` | FILTER: pairwise Cα RMSD + greedy max-min selection | `in`, `chain_field` (default `chain_domain`), `n_select` | `manifest`, `out_dir`, `matrix` |
| `plot_heatmap` | PLOT: render a pairwise RMSD matrix as a heatmap | `matrix` (a saved `rmsd_matrix.csv`) **or** `in` + `chain_field` (compute fresh), `title`, `dpi`, `annotate` | `image` |
| `generate_each` | GENERATE: render one figure per candidate in a manifest | `selection`, `script` (default `pdz_figure.py`), plus any flags to forward | `pngs`, `out_dir` |
| `render` | GENERATE: one-off render for a script that takes a single (input, output) pair instead of looping over a manifest (e.g. `ligand_hotspot_figure.py`) | `script`, `input`, `out`, plus any flags to forward | `image` |
| `render_pair` | GENERATE: one-off render for a script that produces two coupled panels from one shared alignment in a single PyMOL session (e.g. `aligned_pair_figure.py`) instead of one output | `script`, `reference`, `design`, `out_reference`, `out_design`, plus any flags to forward | `reference_image`, `design_image` |
| `render_overlay` | GENERATE: one-off render for a script that takes two structures and produces one combined, superposed panel (e.g. `aligned_overlay_figure.py`) | `script`, `reference`, `design`, `out`, plus any flags to forward | `image` |
| `assemble` | ASSEMBLE: tile images into a montage, optionally rescale | `images`, `rows`, `cols`, `out`, `target_width_in`, `dpi`, `padding`, `bg`, `no_scale` | `image` |
| `assemble_panel_layout` | ASSEMBLE: a wide "problem" panel next to a rows x cols grid of "design" panels — a layout a uniform grid can't express | `left`, `right`, `right_rows`, `right_cols`, `left_width_fraction`, `out`, `target_width_in`, `dpi`, `padding`, `bg`, `no_scale` | `image` |

All worked-example configs and their bash-wrapper counterparts live under
`examples/<name>/`, one subdirectory per example, each runnable from inside
its own directory (see each file's header comment for the exact command):

| Directory | Contents |
|---|---|
| `examples/heatmap/` | `rmsd_heatmap_pipeline.yaml` — the diversity+heatmap demo above |
| `examples/reference_experiment/` | `reference_vs_experiment.yaml` (+ its bash-wrapper twin `run_reference_vs_experiment.sh`) — two independent FIND→FILTER→GENERATE→ASSEMBLE branches (a reference directory, and a production campaign filtered to one best design per target) each reduced to a 5-panel montage, then assembled together into a single comparison figure |
| `examples/diverse_figures/` | `run_diverse_figures.sh` — the bash-wrapper predecessor of the heatmap pipeline (FIND+FILTER+GENERATE+ASSEMBLE over one directory, no YAML) |
| `examples/smallmol/` | `small_molecule_binder_comparison.yaml` — a left "basic problem" panel (the target ligand alone, hot-spot atoms highlighted) next to a 2x2 grid of designs — one best-pLDDT and one median-pLDDT design (`filter_top_n --select top_and_median`) from each of two campaigns (adaptive production vs. nonadaptive reference) — assembled with `assemble_panel_layout`; each design panel carries the problem panel's red/yellow hot-spot coloring through onto its bound peptide via `generate_figure.py --peptide-design-json`. Also the canonical example of the `data:`/`${data.NAME}`/`--root` mechanism — rerun against a different production campaign directory with `--root prod=/path/to/experiment-2`, no YAML edits needed |
| `examples/backbone_fold_overlay/` | `run_backbone_fold_overlay.sh` — bash-wrapper-only example (no YAML/FIND/FILTER, one hand-picked backbone/fold pair): renders a de novo RFDiffusion3 backbone against the downstream AlphaFold2/ColabFold fold prediction it seeds, aligned via `lib/oligomer_align.py`'s best-repeat-window search (the fold predicts the whole oligomeric assembly as one long single chain, not just the one protomer the backbone is) |
| `examples/discontinuous_scaffolds_motif/` | `motif_panels_pipeline.yaml` (+ `resolve_panels.sh`/`run_motif_panels.sh`) — a 1x5 row, one discontinuous-scaffolds design per RESIDUE_ISLAND_COUNT (2-6), each panel built by `motif_superposition_figure.py`: reference ligand (licorice) with the best-passing folded design's protein (cartoon) Kabsch-aligned onto it and its motif hot-spot atoms highlighted (spheres) |
| `examples/discontinuous_scaffolds_motif_single/` | `motif_single_pipeline.yaml` (+ `resolve_panel.sh`/`run_motif_single.sh`) — a single close-up panel for one four-island design: same cartoon+spheres as above, but the ligand as a 50%-transparent surface and the camera rotated so the motif hot spot (not the ligand) faces the viewer, via `motif_superposition_figure.py`'s `--ligand-representation surface` / `--orient-toward hotspot` / `--zoom-target motif` |
| `examples/pdz_design_vs_template/` | `design_vs_template_pipeline.yaml` (+ `resolve_template.py`/`resolve_template.sh`/`run_design_vs_template.sh`) — two views of the pdzbinder production campaign's single highest-confidence design vs. its origin template, both built on the same C-terminal peptide-backbone Kabsch fit (`lib/peptide_align.py`): a two-panel row (`aligned_pair_figure.py`, domain cyan / peptide green in both panels) and a single superposed overlay (`aligned_overlay_figure.py`, template domain cyan / template peptide green / design domain red / design peptide yellow) |
| `examples/discontinuous_scaffolds_anchor_progression/` | `anchor_progression_pipeline.yaml` (+ `resolve_progression.sh`/`run_anchor_progression.sh`) — a 1x2 row for one input motif (M0157_1qh5) showing its design progression across the one pipeline run (of several) whose adaptive fold-redesign lineage ends in a passing model: left panel is the initial (root) design with its anchor residues highlighted, right panel is the final (passing redesign) design with the same residues highlighted to show the portion of the motif that was subject to anchoring — both anchor and non-anchor hot-spot residues shown as colored cartoon segments (burnt orange / purple), no spheres — via `motif_superposition_figure.py`'s new `anchor_positions`/`--anchor-color`/`--motif-representation cartoon` support and `lib/rfd3_motif_select.py`'s new `anchor_chai_positions`. `anchor_progression_pipeline_M0349_1e3v.yaml` (+ `resolve_progression.sh M0349_1e3v`/`run_anchor_progression_M0349_1e3v.sh`) is a second worked model, picked as the one with the most root-generation anchor residues (2, vs. M0157_1qh5's 1) among every model whose redesign lineage ends in a passing model |

## Function reference

| Script | Role |
|---|---|
| `find_structures_flat.py` | FIND — glob a flat directory of PDBs into a manifest |
| `find_structures_campaign.py` | FIND — walk a nested production-campaign layout, attach Boltz confidence scores, normalize chain ids |
| `find_structures_smallmol.py` | FIND — discover completed small-molecule-binder designs from an ImpressBasePipeline campaign, pairing each AF2 fold-prediction with its preceding FastRelax protein+ligand complex |
| `filter_best_score.py` | FILTER — one best-scoring candidate per group key |
| `filter_top_n.py` | FILTER — top N candidates overall by a score field, no grouping; or (`--select top_and_median`) a fixed best-scorer-plus-median-scorer pick |
| `filter_diversity.py` | FILTER — pairwise-dissimilar subset by chain-CA RMSD (wraps `lib/rmsd.py`) |
| `plot_rmsd_heatmap.py` | PLOT — render a pairwise RMSD matrix as a heatmap PNG |
| `pdz_figure.py` | GENERATE — auto-oriented PDZ-domain/peptide complex figure (fixed chain A/B, pink/lime-green) |
| `generate_figure.py` | GENERATE — the same auto-oriented figure, generalized to any chain ids/colors; optionally colors the peptide chain by an RFDiffusion3 binder-design spec's hot-spot split (`--peptide-design-json`) instead of one flat color |
| `ligand_hotspot_figure.py` | GENERATE — render a target ligand alone, split into hot-spot (red) vs. rest (yellow) atoms per an RFDiffusion3 binder-design spec |
| `backbone_fold_overlay_figure.py` | GENERATE — render a de novo RFDiffusion3 backbone (one protomer) against the downstream AlphaFold2/ColabFold fold prediction of the whole oligomeric assembly it seeds, aligned via `lib/oligomer_align.py`'s best-repeat-window search; solo views of each, two aligned overlap panels (representations/colors swapped between passes), plus optional fold-alone ribbon/secondary-structure-colored extras |
| `motif_superposition_figure.py` | GENERATE — render a folded design's protein (cartoon) and its RFD3 motif hot-spot atoms (spheres by default, or a colored cartoon segment via `--motif-representation cartoon`), Kabsch-aligned onto a reference active-site structure whose ligand is drawn as licorice (default) or a transparent surface; camera can face the ligand (default) or the hot spot itself; an optional per-panel `anchor_positions` list always renders that subset as spheres in a second color (default burnt orange), regardless of `--motif-representation` |
| `aligned_pair_figure.py` | GENERATE — render a design-template ("reference") complex and the resulting ("design") complex as two separate cartoon panels (domain + peptide, one color each) sharing a coordinate frame, the reference's peptide Kabsch-fit onto the design's peptide's C-terminal residues |
| `aligned_overlay_figure.py` | GENERATE — the same reference/design C-terminal peptide-backbone alignment as `aligned_pair_figure.py`, rendered as one combined panel with both structures superposed (four independent cartoon colors: reference domain/peptide, design domain/peptide) |
| `anchor_progression_oriented_figure.py` | GENERATE — two-panel root-vs-redesign anchor-progression figure sharing one fixed camera orientation (computed once from the root design's motif plane, then reapplied verbatim to the redesign panel) instead of `motif_superposition_figure.py`'s per-panel independent orientation; flat protein color plus colored cartoon/sphere segments for hot-spot/anchor residues |
| `anchor_progression_zoned_figure.py` | GENERATE — the same shared-orientation two-panel anchor-progression figure as `anchor_progression_oriented_figure.py`, but coloring the protein by distance from the motif (near = opaque cartoon, far = translucent surface) instead of one flat color |
| `montage_figures.py` | ASSEMBLE — tile images into a grid, optionally rescale to a target width |
| `assemble_panel_layout.py` | ASSEMBLE — a wide left panel next to a rows x cols grid of right-hand panels, a layout `montage_figures.py`'s uniform grid can't express |
| `run_pipeline.py` | Generic runner that dispatches a YAML pipeline's steps to the scripts above |
| `pdz_pairwise_rmsd.py` | Standalone predecessor of `filter_diversity.py` (directory glob instead of a manifest, hard-coded chain A) — kept for reference |
| `rmsd_heatmap_notebook.py` | Interactive marimo notebook: documents the YAML format and runs find → filter_diversity → plot_rmsd_heatmap live |
| `lib/manifest.py` | The candidate-manifest read/write contract shared by every step |
| `lib/rmsd.py` | PyMOL-backed chain-CA RMSD core: load structures, compute the pairwise matrix, greedy max-min selection, chain-id-agreement check, matrix CSV I/O |
| `lib/kabsch.py` | Point-correspondence rigid-body superposition (Kabsch algorithm) for two paired Nx3 coordinate arrays |
| `lib/pdb_normalize.py` | Detect and rewrite non-standard multi-character PDB chain ids |
| `lib/ligand_select.py` | Atom-name-based ligand selection from an RFDiffusion3 binder-design spec, for splitting one hetero-residue into sub-groups |
| `lib/rfd3_motif_select.py` | RFD3 structured design-spec selection: contig -> chai-position mapping, protein-motif-residue vs. ligand atom-name split, and analysis.py `anchor_residues` CSV parsing -> anchor chai-position sets |
| `lib/peptide_align.py` | C-terminal peptide-backbone Kabsch alignment: match one chain's last N residues to another's by position (N -> C order), independent of residue identity or an external contig map |
| `lib/oligomer_align.py` | cealign a single protomer onto its best-matching repeat window within a downstream single-chain oligomeric fold prediction, sliding a CA-count-sized window and keeping whichever scores best |
| `lib/design_spec.py` | Atom-name-CSV parsing and spec-relative-path resolution shared by `lib/ligand_select.py` and `lib/rfd3_motif_select.py` |
| `lib/imgtrim.py` | Crop a rendered figure to its content bounding box before tiling into a montage |
| `lib/montage.py` | Grid-tiling (`build_montage`) and target-width scaling (`scale_to_width`) shared by `montage_figures.py` and `assemble_panel_layout.py` |
| `lib/geometry.py` | Pure-numpy orientation math: Rodrigues rotation (`rotation_to_align`) and PCA long-axis/plane fitting, no PyMOL dependency |
| `lib/pymol_scene.py` | PyMOL session bootstrap, the shared "AOShiny" material settings, ray-trace-and-save, and the `--width/--height/--dpi/--bg` flag helper — used by every GENERATE script |
| `lib/orient.py` | The actual scene-orientation routines (long-axis-vertical with optional twist; look-down-on-plane) built on `lib/geometry.py` + `lib/kabsch.py`, used by all four GENERATE scripts |

## Requirements

A Python 3 environment with PyMOL's Python API importable (`import pymol`),
plus `numpy`, `PyYAML`, `Pillow`, and `matplotlib`. `pytest` is only needed
to run the test suite (`pytest tests/`), `marimo` only to run the notebook.
There's no dependency manifest — these must already be on the interpreter's
`PATH`/`site-packages`.
