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
python3 run_pipeline.py rmsd_heatmap_pipeline.yaml
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
| `filter_best_score` | FILTER: keep the best-scoring candidate per group (e.g. one design per target) | `in`, `key`, `score_field`, `mode` (`max`\|`min`) | `manifest` |
| `filter_diversity` | FILTER: pairwise Cα RMSD + greedy max-min selection | `in`, `chain_field` (default `chain_domain`), `n_select` | `manifest`, `out_dir`, `matrix` |
| `plot_heatmap` | PLOT: render a pairwise RMSD matrix as a heatmap | `matrix` (a saved `rmsd_matrix.csv`) **or** `in` + `chain_field` (compute fresh), `title`, `dpi`, `annotate` | `image` |
| `generate_each` | GENERATE: render one figure per candidate in a manifest | `selection`, `script` (default `pdz_figure.py`), plus any flags to forward | `pngs`, `out_dir` |
| `assemble` | ASSEMBLE: tile images into a montage, optionally rescale | `images`, `rows`, `cols`, `out`, `target_width_in`, `dpi`, `padding`, `bg`, `no_scale` | `image` |

`reference_vs_experiment.yaml` is a larger worked example: two independent
SELECT→GENERATE→ASSEMBLE branches (a reference directory, and a production
campaign filtered to one best design per target) each reduced to a 5-panel
montage, then assembled together into a single comparison figure.

## Function reference

| Script | Role |
|---|---|
| `find_structures_flat.py` | FIND — glob a flat directory of PDBs into a manifest |
| `find_structures_campaign.py` | FIND — walk a nested production-campaign layout, attach Boltz confidence scores, normalize chain ids |
| `filter_best_score.py` | FILTER — one best-scoring candidate per group key |
| `filter_diversity.py` | FILTER — pairwise-dissimilar subset by chain-CA RMSD (wraps `lib/rmsd.py`) |
| `plot_rmsd_heatmap.py` | PLOT — render a pairwise RMSD matrix as a heatmap PNG |
| `pdz_figure.py` | GENERATE — auto-oriented PDZ-domain/peptide complex figure (fixed chain A/B, pink/lime-green) |
| `generate_figure.py` | GENERATE — the same auto-oriented figure, generalized to any chain ids/colors |
| `montage_figures.py` | ASSEMBLE — tile images into a grid, optionally rescale to a target width |
| `run_pipeline.py` | Generic runner that dispatches a YAML pipeline's steps to the scripts above |
| `pdz_pairwise_rmsd.py` | Standalone predecessor of `filter_diversity.py` (directory glob instead of a manifest, hard-coded chain A) — kept for reference |
| `rmsd_heatmap_notebook.py` | Interactive marimo notebook: documents the YAML format and runs find → filter_diversity → plot_rmsd_heatmap live |
| `lib/manifest.py` | The candidate-manifest read/write contract shared by every step |
| `lib/rmsd.py` | PyMOL-backed chain-CA RMSD core: load structures, compute the pairwise matrix, greedy max-min selection |
| `lib/pdb_normalize.py` | Detect and rewrite non-standard multi-character PDB chain ids |

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

# Or run a whole pipeline
python3 run_pipeline.py rmsd_heatmap_pipeline.yaml

# Or explore interactively
marimo edit rmsd_heatmap_notebook.py
```
