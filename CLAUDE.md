# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

molvizgen is a small collection of CLI scripts for building a
SELECT → GENERATE → ASSEMBLE pipeline over sets of protein structure (PDB)
files: enumerate candidate structures, filter them down (by score or by
structural diversity), render publication-style PyMOL figures or an RMSD
heatmap, and assemble the results into montages. There is no test suite,
linter, or build step — each script is a standalone `argparse` CLI, runnable
directly.

This repo holds only code. It operates on PDB files that live elsewhere
(pass a directory via `--dir`/`dir:`); it has no data of its own.

## Commands

Everything requires a Python environment with the PyMOL Python API
importable (`import pymol`), plus `numpy`, `PyYAML`, `Pillow`, `matplotlib`,
and `marimo` (only for the notebook). There's no dependency manifest in this
repo — these must already be present in whatever `python3` is on `PATH`.

Run a single step directly, e.g.:
```
python3 find_structures_flat.py --dir /path/to/pdbs --out candidates.json
python3 filter_diversity.py --in candidates.json --out-dir analysis
python3 plot_rmsd_heatmap.py --matrix analysis/rmsd_matrix.csv --out heatmap.png
```

Run a full YAML-described pipeline (see `run_pipeline.py`'s docstring and
`examples/heatmap/rmsd_heatmap_pipeline.yaml` /
`examples/reference_experiment/reference_vs_experiment.yaml` /
`examples/smallmol/small_molecule_binder_comparison.yaml` for worked
examples). Each example config lives in its own `examples/<name>/`
subdirectory and is meant to be run from inside it, e.g.:
```
cd examples/heatmap && python3 ../../run_pipeline.py rmsd_heatmap_pipeline.yaml
```

Open the interactive notebook (documents the YAML format and runs
find → filter_diversity → plot_rmsd_heatmap live):
```
marimo edit rmsd_heatmap_notebook.py
```

There are also two hand-written bash wrappers,
`examples/diverse_figures/run_diverse_figures.sh` and
`examples/reference_experiment/run_reference_vs_experiment.sh`, that chain
the same underlying scripts without going through `run_pipeline.py`/YAML —
useful as a second reference for how the scripts compose, or for one-off
runs. Both locate the repo's Python scripts via a `REPO_ROOT` computed as
two directories up from their own location (`examples/<name>/../..`), since
the scripts themselves stay at the repo root.

## Architecture

**The manifest is the contract.** Every step in the pipeline (`lib/manifest.py`)
passes a list of candidate dicts as JSON: `{"candidates": [{"id", "pdb_path",
"chain_domain", "chain_peptide", ...arbitrary scoring fields}, ...]}`.
`pdb_path` is always absolute, so no step ever has to guess a file location
from an id + directory convention. FIND steps produce a manifest from
scratch; FILTER steps consume one manifest and produce a smaller/reordered
one; GENERATE steps consume a manifest and produce images.

**`run_pipeline.py` is a generic dispatcher, not a workflow engine.** A YAML
pipeline is `out_dir` + an ordered list of `{name, kind, args}` steps. Each
step's `kind` maps to a handler function (`HANDLERS` dict at the bottom of
the file) that shells out to one of the CLI scripts in this directory,
translating `args` keys 1:1 into `--flag` arguments
(`flags_from_args`/`flag`). A later step references an earlier one's declared
output with `${step_name.field}` (`resolve()` / `REF_PATTERN`) — this only
matches when the *entire* string is the placeholder; it cannot be embedded
inside a larger string, so a step must expose exactly the field a later step
needs (e.g. `filter_diversity` declares `manifest`, `out_dir`, *and* `matrix`
so a downstream `plot_heatmap` step can reference the CSV path directly).
Adding a new step kind means writing one `handle_<kind>` function and
registering it in `HANDLERS`.

**`lib/rmsd.py` is the structural-diversity core**, used by both
`filter_diversity.py` (the manifest-based, general FILTER step) and the
older `pdz_pairwise_rmsd.py` (a self-contained, directory-glob-based
predecessor kept for reference). It loads chain-CA-only PyMOL objects
(`load_ca_objects`), computes the full symmetric pairwise RMSD matrix via
`cealign` (`pairwise_rmsd_matrix`), and picks an `n`-sized
maximally-dissimilar subset via greedy farthest-point search
(`greedy_max_min_selection`: seed with the single most dissimilar pair, then
repeatedly add whichever remaining candidate has the largest *minimum*
distance to the already-selected set). `plot_rmsd_heatmap.py` renders that
matrix (read from a saved CSV, or computed fresh from a manifest) as a PNG.

**Chain-id normalization is a real gotcha handled in `lib/pdb_normalize.py`.**
The strict PDB format reserves one fixed column for chain ID; some upstream
pipelines (observed: Boltz output) instead write multi-character,
whitespace-separated tokens there (e.g. `"pdz"`/`"pep"`). Fixed-column
readers — including PyMOL — silently truncate both to `"p"`, colliding
without erroring. `find_structures_campaign.py` detects this
(`detect_chain_tokens`) and, when found, rewrites affected files into cached,
standards-compliant single-character-chain copies (`normalize_pdb_chains`)
before they ever reach a FILTER/GENERATE step — downstream code always deals
in plain `chain A`/`chain B` selectors.

**Figure rendering (`pdz_figure.py`, its generalized twin
`generate_figure.py`) auto-orients structures** rather than using a fixed
camera: it PCA's the domain chain's CA atoms to find its long axis, rotates
that to vertical, and picks the remaining twist angle so the peptide chain
faces the camera — so every rendered figure is comparably framed regardless
of the structure's original coordinate frame. `pdz_figure.py` hard-codes the
PDZ/peptide chain-A/chain-B convention and pink/lime-green coloring;
`generate_figure.py` is the same logic with chain ids and colors as flags,
for structures that don't use that convention. `ligand_hotspot_figure.py`
applies the same long-axis-to-vertical orientation to a lone ligand (no
second chain to twist toward the camera, since there's no protein in that
figure at all).

**Atom-name selection (`lib/ligand_select.py`) is a second selection
convention alongside chain ids**, for figures that need to split a single
hetero-residue into sub-groups rather than select whole chains — e.g. the
"hot spot" atoms an RFDiffusion3 binder-design spec (an ImpressBasePipeline
`*_binder_design.json` file) designates as buried/targeted vs. the rest of
the ligand. `load_binder_design_spec` reads that JSON and resolves its
`input` PDB path relative to the JSON's own directory (same convention every
pipeline here uses); `atom_name_selection` builds the PyMOL `name`-based
selection expression from the resulting atom-name lists.
`ligand_hotspot_figure.py` is the GENERATE step built on top of it.

**`find_structures_smallmol.py` mirrors an external pipeline's layout and
state machine** (ImpressBasePipeline's `SmallMoleculeBindingPipeline` — see
that project's own CLAUDE.md for the authoritative state machine). A
campaign root holds flat `p<N>/` directories, each a numbered sequence of
`<taskcount>_<taskname>/{in,out}/` task directories. A "completed design" is
one that reached an `alphafold` task with a computable mean pLDDT; because
AF2/ColabFold doesn't co-fold small molecules, the AF2 output itself is
protein-only, so this step walks back to the `fastrelax` task that
immediately precedes it (observed convention: consecutive task numbers) to
find the actual protein+ligand complex PDB worth rendering.

`montage_figures.py` tiles a list of images into an evenly-spaced grid and
optionally rescales the whole montage to a target width — intermediate
montages in a larger assembly should pass `--no-scale` and only the final
assemble step should scale, so repeated resizing doesn't degrade image
quality. Before tiling, both it and `assemble_panel_layout.py` trim each
source image to its content bounding box via `lib/imgtrim.py`
(`trim_to_content`, `--no-trim` to skip) — every GENERATE script rays-traces
onto a fixed canvas with its own zoom buffer, so raw PNGs otherwise carry
inconsistent background margin that's most visible when panels have
different aspect ratios. `assemble_panel_layout.py` handles a layout a
uniform `rows x cols` grid can't express: a wide "problem" panel spanning
the full height of an adjacent `rows x cols` grid of "design" panels (built
by calling `montage_figures.py`'s `build_montage()` directly), the left
panel's width computed as a fraction of the *final* canvas.
