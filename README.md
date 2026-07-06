# molvizgen

Command-line tools for screening a set of protein structure (PDB) files down
to a representative subset and turning them into figures: a
**SELECT → GENERATE → ASSEMBLE** pipeline built from small, composable
scripts, driven either directly or by a YAML pipeline description.

molvizgen is a *tool* repo, not a *data* repo — every script takes a
directory (or a manifest) of PDB files as input via a flag; there are no
structures checked in here.

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
lone ligand, since there's no second chain to twist toward the camera.

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
small pad) before tiling; both `montage_figures.py` and
`assemble_panel_layout.py` do this by default (`--no-trim` to skip).

## The pipeline

A pipeline is a YAML file: a base output directory plus an ordered list of
named steps, each with a `kind` and `args`:

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
| `assemble` | ASSEMBLE: tile images into a montage, optionally rescale | `images`, `rows`, `cols`, `out`, `target_width_in`, `dpi`, `padding`, `bg`, `no_scale` | `image` |
| `assemble_panel_layout` | ASSEMBLE: a wide "problem" panel next to a rows x cols grid of "design" panels — a layout a uniform grid can't express | `left`, `right`, `right_rows`, `right_cols`, `left_width_fraction`, `out`, `target_width_in`, `dpi`, `padding`, `bg`, `no_scale` | `image` |

All worked-example configs and their bash-wrapper counterparts live under
`examples/<name>/`, one subdirectory per example, each runnable from inside
its own directory (see each file's header comment for the exact command):

| Directory | Contents |
|---|---|
| `examples/heatmap/` | `rmsd_heatmap_pipeline.yaml` — the diversity+heatmap demo above |
| `examples/reference_experiment/` | `reference_vs_experiment.yaml` (+ its bash-wrapper twin `run_reference_vs_experiment.sh`) — two independent SELECT→GENERATE→ASSEMBLE branches (a reference directory, and a production campaign filtered to one best design per target) each reduced to a 5-panel montage, then assembled together into a single comparison figure |
| `examples/diverse_figures/` | `run_diverse_figures.sh` — the bash-wrapper predecessor of the heatmap pipeline (FIND+FILTER+GENERATE+ASSEMBLE over one directory, no YAML) |
| `examples/smallmol/` | `small_molecule_binder_comparison.yaml` — a left "basic problem" panel (the target ligand alone, hot-spot atoms highlighted) next to a 2x2 grid of top-pLDDT designs drawn from two campaigns (adaptive production vs. nonadaptive reference), assembled with `assemble_panel_layout` |

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
| `montage_figures.py` | ASSEMBLE — tile images into a grid, optionally rescale to a target width |
| `assemble_panel_layout.py` | ASSEMBLE — a wide left panel next to a rows x cols grid of right-hand panels, a layout `montage_figures.py`'s uniform grid can't express |
| `run_pipeline.py` | Generic runner that dispatches a YAML pipeline's steps to the scripts above |
| `pdz_pairwise_rmsd.py` | Standalone predecessor of `filter_diversity.py` (directory glob instead of a manifest, hard-coded chain A) — kept for reference |
| `rmsd_heatmap_notebook.py` | Interactive marimo notebook: documents the YAML format and runs find → filter_diversity → plot_rmsd_heatmap live |
| `lib/manifest.py` | The candidate-manifest read/write contract shared by every step |
| `lib/rmsd.py` | PyMOL-backed chain-CA RMSD core: load structures, compute the pairwise matrix, greedy max-min selection |
| `lib/pdb_normalize.py` | Detect and rewrite non-standard multi-character PDB chain ids |
| `lib/ligand_select.py` | Atom-name-based ligand selection from an RFDiffusion3 binder-design spec, for splitting one hetero-residue into sub-groups |
| `lib/imgtrim.py` | Crop a rendered figure to its content bounding box before tiling into a montage |

## Requirements

A Python 3 environment with PyMOL's Python API importable (`import pymol`),
plus `numpy`, `PyYAML`, `Pillow`, and `matplotlib`. `marimo` is only needed
to run the notebook. There's no dependency manifest — these must already be
on the interpreter's `PATH`/`site-packages`.

## Usage

```bash
# Run one step directly
python3 find_structures_flat.py --dir /path/to/pdbs --out candidates.json
python3 filter_diversity.py --in candidates.json --out-dir analysis
python3 plot_rmsd_heatmap.py --matrix analysis/rmsd_matrix.csv --out heatmap.png

# Or run a whole pipeline (see examples/<name>/ for worked configs)
cd examples/heatmap && python3 ../../run_pipeline.py rmsd_heatmap_pipeline.yaml

# Or explore interactively
marimo edit rmsd_heatmap_notebook.py
```
