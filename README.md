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

**Pairwise RMSD diversity.** Given a set of structures that are all "the same
kind of thing" (e.g. the same domain solved/predicted many times), the
useful question is usually not "what's the single best one" but "give me a
handful that are maximally *different* from each other" — a diverse panel
for a figure, or a diverse benchmark set. molvizgen answers this by:
1. aligning every pair of structures on a chosen chain's Cα atoms (PyMOL
   `cealign`) and recording the RMSD, producing a full symmetric matrix;
2. greedily selecting a subset that maximizes the *minimum* pairwise RMSD
   within it (farthest-point/max-min selection) — seed with the single most
   dissimilar pair, then repeatedly add whichever remaining structure is
   farthest from everything already picked.

**Chain-id normalization.** The PDB format reserves a single fixed column for
chain ID. Some upstream pipelines (observed: Boltz) instead write
multi-character tokens there as whitespace-separated fields (`"pdz"`/`"pep"`).
Naive fixed-column readers — PyMOL included — silently truncate both to the
same single letter, corrupting chain-based selection without raising an
error. molvizgen detects this up front and rewrites affected files into
standards-compliant single-character-chain copies before anything downstream
touches them.

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

**Atom-name selection.** Chain ids select whole chains; some figures need to
split a *single* hetero-residue into sub-groups instead (e.g. the hot-spot
atoms an RFDiffusion3 binder-design spec targets vs. the rest of the
ligand). `lib/ligand_select.py` reads that split from an
`*_binder_design.json` spec and builds a PyMOL `name`-based selection
expression, complementing the chain-id convention used elsewhere.

**Content-trimmed montages.** Every GENERATE script ray-traces onto a
fixed-size canvas with its own zoom buffer, so raw PNGs carry inconsistent
background margin — most visible when panels have different aspect ratios
(e.g. a tall thin ligand next to a roughly square protein cartoon).
`lib/imgtrim.py` crops each source image to its content bounding box (plus a
small pad) before tiling; `lib/montage.py`'s `build_montage` does this by
default (`--no-trim` to skip), shared by both `montage_figures.py` and
`assemble_panel_layout.py` rather than each carrying its own tiling logic.

**Point-correspondence alignment (Kabsch).** Some figures already know an
explicit atom-to-atom correspondence between two structures — e.g. an
RFDiffusion3 contig says reference residue (chain, resnum) is the same
residue as sequence position N in a folded design — rather than needing
`cealign`'s own sequence/structure alignment to find one. `lib/kabsch.py`
computes the least-squares rotation/translation for that case, used by
`motif_superposition_figure.py` to place a folded design's motif into a
reference structure's coordinate frame.

**RFD3 design-spec selection.** A third selection convention, alongside
chain ids and the buried/exposed ligand atom-name split: `lib/rfd3_motif_select.py`
reads the per-model design spec used by the discontinuous-scaffolds
RFDiffusion3 benchmark (`{model: {input, ligand, contig, select_fixed_atoms}}`),
splitting `select_fixed_atoms` into protein motif residues vs. ligand atoms
and mapping each motif residue to its 1-indexed position in any fold
generated from that contig (including redesign generations — see the
module docstring).

**C-terminal peptide-backbone alignment.** A fourth selection/alignment
convention: `lib/peptide_align.py` Kabsch-superposes one peptide chain's
last N residues onto another's, matched purely by position (N -> C order)
rather than residue identity or an external contig map — the right tool
when a short peptide fragment (e.g. a PDZ domain's conserved C-terminal
binding motif, used as a design template) and a longer peptide produced
downstream share that same C-terminal motif but otherwise differ freely in
length or sequence. `lib/peptide_align.py`'s `load_and_align_pair` bundles
the shared recipe both figure scripts built on it need: load a `design`
complex oriented with the standard long-axis-vertical convention, and a
`reference` complex whose coordinates are rigidly transformed so its
peptide's C-terminal motif lands on the design's. `aligned_pair_figure.py`
renders that pair as **two** separate, still-comparable panels (meant to be
tiled side by side); `aligned_overlay_figure.py` renders the identical
alignment as **one** combined, superposed panel instead — closer to
`motif_superposition_figure.py`'s single-panel convention, but comparing two
full complexes (four independent colors) rather than a minimal reference
extract's ligand plus a design's motif spheres.

## The pipeline

A pipeline is a YAML file: a base output directory plus an ordered list of
named steps, each with a `kind` and `args`. Every `kind` belongs to one of
five stages: **FIND** (enumerate a directory or campaign layout into a
manifest), **FILTER** (reduce a manifest by score or structural diversity),
**PLOT** (render a saved RMSD matrix as a heatmap), **GENERATE** (render
PyMOL figures, per-candidate or one-off), and **ASSEMBLE** (tile rendered
figures into a montage or panel layout) — see [Step kinds](#step-kinds) below
for the full list. The example below chains FIND → FILTER → PLOT:

```yaml
out_dir: analysis/rmsd_heatmap_demo

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
subdirectory, `out_dir/<step name>/`.

### Step kinds

| `kind` | Purpose | Key `args` | Declares |
|---|---|---|---|
| `find_flat` | FIND: enumerate every PDB file in a directory | `dir`, `glob` (default `*.pdb`), `chain_domain`, `chain_peptide` | `manifest` |
| `find_campaign` | FIND: recursively discover structures from a production campaign layout, attach confidence scores, normalize non-standard chain ids | `prod_root`, `groups`, `raw_chain_domain`, `raw_chain_peptide`, `normalized_dir` | `manifest` |
| `find_smallmol` | FIND: discover completed small-molecule-binder designs from an ImpressBasePipeline campaign layout, pairing each AF2 fold-prediction with its preceding FastRelax protein+ligand complex | `campaign_root`, `campaign_label`, `chain_domain`, `chain_peptide` | `manifest` |
| `filter_best_score` | FILTER: keep the best-scoring candidate per group (e.g. one design per target) | `in`, `key`, `score_field`, `mode` (`max`\|`min`) | `manifest` |
| `filter_top_n` | FILTER: keep the top N candidates overall by a score field, no grouping | `in`, `score_field`, `n`, `mode` (`max`\|`min`) | `manifest` |
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
| `examples/smallmol/` | `small_molecule_binder_comparison.yaml` — a left "basic problem" panel (the target ligand alone, hot-spot atoms highlighted) next to a 2x2 grid of top-pLDDT designs drawn from two campaigns (adaptive production vs. nonadaptive reference), assembled with `assemble_panel_layout` |
| `examples/discontinuous_scaffolds_motif/` | `motif_panels_pipeline.yaml` (+ `resolve_panels.sh`/`run_motif_panels.sh`) — a 1x5 row, one discontinuous-scaffolds design per RESIDUE_ISLAND_COUNT (2-6), each panel built by `motif_superposition_figure.py`: reference ligand (licorice) with the best-passing folded design's protein (cartoon) Kabsch-aligned onto it and its motif hot-spot atoms highlighted (spheres) |
| `examples/discontinuous_scaffolds_motif_single/` | `motif_single_pipeline.yaml` (+ `resolve_panel.sh`/`run_motif_single.sh`) — a single close-up panel for one four-island design: same cartoon+spheres as above, but the ligand as a 50%-transparent surface and the camera rotated so the motif hot spot (not the ligand) faces the viewer, via `motif_superposition_figure.py`'s `--ligand-representation surface` / `--orient-toward hotspot` / `--zoom-target motif` |
| `examples/pdz_design_vs_template/` | `design_vs_template_pipeline.yaml` (+ `resolve_template.py`/`resolve_template.sh`/`run_design_vs_template.sh`) — two views of the pdzbinder production campaign's single highest-confidence design vs. its origin template, both built on the same C-terminal peptide-backbone Kabsch fit (`lib/peptide_align.py`): a two-panel row (`aligned_pair_figure.py`, domain cyan / peptide green in both panels) and a single superposed overlay (`aligned_overlay_figure.py`, template domain cyan / template peptide green / design domain red / design peptide purple) |

## Function reference

| Script | Role |
|---|---|
| `find_structures_flat.py` | FIND — glob a flat directory of PDBs into a manifest |
| `find_structures_campaign.py` | FIND — walk a nested production-campaign layout, attach Boltz confidence scores, normalize chain ids |
| `find_structures_smallmol.py` | FIND — discover completed small-molecule-binder designs from an ImpressBasePipeline campaign, pairing each AF2 fold-prediction with its preceding FastRelax protein+ligand complex |
| `filter_best_score.py` | FILTER — one best-scoring candidate per group key |
| `filter_top_n.py` | FILTER — top N candidates overall by a score field, no grouping |
| `filter_diversity.py` | FILTER — pairwise-dissimilar subset by chain-CA RMSD (wraps `lib/rmsd.py`) |
| `plot_rmsd_heatmap.py` | PLOT — render a pairwise RMSD matrix as a heatmap PNG |
| `pdz_figure.py` | GENERATE — auto-oriented PDZ-domain/peptide complex figure (fixed chain A/B, pink/lime-green) |
| `generate_figure.py` | GENERATE — the same auto-oriented figure, generalized to any chain ids/colors |
| `ligand_hotspot_figure.py` | GENERATE — render a target ligand alone, split into hot-spot (red) vs. rest (yellow) atoms per an RFDiffusion3 binder-design spec |
| `motif_superposition_figure.py` | GENERATE — render a folded design's protein (cartoon) and its RFD3 motif hot-spot atoms (spheres), Kabsch-aligned onto a reference active-site structure whose ligand is drawn as licorice (default) or a transparent surface; camera can face the ligand (default) or the hot spot itself |
| `aligned_pair_figure.py` | GENERATE — render a design-template ("reference") complex and the resulting ("design") complex as two separate cartoon panels (domain + peptide, one color each) sharing a coordinate frame, the reference's peptide Kabsch-fit onto the design's peptide's C-terminal residues |
| `aligned_overlay_figure.py` | GENERATE — the same reference/design C-terminal peptide-backbone alignment as `aligned_pair_figure.py`, rendered as one combined panel with both structures superposed (four independent cartoon colors: reference domain/peptide, design domain/peptide) |
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
| `lib/rfd3_motif_select.py` | RFD3 structured design-spec selection: contig -> chai-position mapping, protein-motif-residue vs. ligand atom-name split |
| `lib/peptide_align.py` | C-terminal peptide-backbone Kabsch alignment: match one chain's last N residues to another's by position (N -> C order), independent of residue identity or an external contig map |
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
