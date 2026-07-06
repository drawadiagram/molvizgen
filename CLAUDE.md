# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

molvizgen is a small collection of CLI scripts for building a
SELECT → GENERATE → ASSEMBLE pipeline over sets of protein structure (PDB)
files: enumerate candidate structures, filter them down (by score or by
structural diversity), render publication-style PyMOL figures or an RMSD
heatmap, and assemble the results into montages. There is no linter or build
step — each script is a standalone `argparse` CLI, runnable directly. A
`pytest` suite under `tests/` covers `lib/`'s shared logic and a CLI smoke
test per pipeline stage (see Commands below).

This repo holds only code. It operates on PDB files that live elsewhere
(pass a directory via `--dir`/`dir:`); it has no data of its own.

## Commands

Everything requires a Python environment with the PyMOL Python API
importable (`import pymol`), plus `numpy`, `PyYAML`, `Pillow`, `matplotlib`,
`pytest` (only for the test suite), and `marimo` (only for the notebook).
There's no dependency manifest in this repo — these must already be present
in whatever `python3` is on `PATH`.

Run the test suite:
```
pytest tests/
```
Tests run against small fixture PDBs/manifests under `tests/fixtures/`
(including a real PDZ-domain/peptide complex and an unrelated fold, both
extracted from local reference structures) rather than the large external
directories the `examples/` pipelines point at, so the suite never depends
on data outside the repo.

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

**That auto-orientation logic, and the PyMOL rendering boilerplate every
GENERATE script needs, live in `lib/`, not in each figure script.**
`lib/geometry.py` is the pure-numpy half (no PyMOL calls): `rotation_to_align`
(Rodrigues' formula) and the PCA helpers `long_axis_pca`/`plane_pca`.
`lib/pymol_scene.py` is the PyMOL-coupled half: importing it launches the
headless PyMOL session as a side effect (the same "import this instead of
calling `pymol.finish_launching` yourself" idiom `lib/rmsd.py` also uses),
plus `apply_material_aoshiny` (the shiny/AO/no-shadows settings),
`ray_trace_and_save`, and `add_render_flags` for the `--width/--height/--dpi/--bg`
flags every figure script exposes identically. `lib/orient.py` builds the
actual scene-orientation routines on top of both: `orient_long_axis_vertical`
is the one function behind all three "long axis vertical, optionally twist
something else toward the camera" conventions above (`face_sel=None` skips
the twist, which is what gives `ligand_hotspot_figure.py`'s no-twist
behavior); `orient_look_down_on_plane` (used by
`motif_superposition_figure.py`'s `--orient-toward hotspot`, see below) is a
distinct plane-normal-to-camera algorithm, not a variant of the long-axis
one. Every GENERATE script (`pdz_figure.py`, `generate_figure.py`,
`ligand_hotspot_figure.py`, `motif_superposition_figure.py`) is a thin
domain-specific layer over these three modules — new figure scripts should
extend `lib/orient.py`'s functions with new parameters rather than
re-deriving this math locally.

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

**`lib/rfd3_motif_select.py` is a third selection convention**, for the
per-model design spec used by the discontinuous-scaffolds RFDiffusion3
benchmark (`{model: {input, ligand, contig, select_fixed_atoms}}` —
see that project's `scripts/mcsa_41-N.json`). `select_fixed_atoms` mixes
protein motif-residue keys ("A1054" = chain + PDB resnum) and ligand
atom-name-split keys (a resn matching a token in "ligand") in one dict;
`split_fixed_atoms` separates them. `contig_to_chai_positions` walks the
RFD3 contig string to map each motif residue to its 1-indexed sequence
position in any fold generated from that contig — including a *redesign*
descendant of it, since a redesign's contig relabels which residues anchor
the scaffold but preserves every free-residue run-length, so the position
arithmetic for a given original motif residue comes out the same regardless
of how many redesign generations produced the fold being examined (see the
discontinuous-scaffolds project's own CLAUDE.md, "Adaptive branching for
failing models", for that state machine).

**`lib/kabsch.py` complements `lib/rmsd.py`'s cealign-based alignment with
point-correspondence alignment**: given two paired Nx3 coordinate arrays
(no sequence/structure alignment step needed, because the correspondence is
already known — e.g. via `contig_to_chai_positions`), `kabsch_fit` returns
the least-squares rotation/translation. `motif_superposition_figure.py`
(GENERATE) composes these two new modules: it Kabsch-aligns a folded
design's motif backbone atoms onto a reference active-site structure (an
M-CSA-style minimal extract — just the fixed motif residues and ligand, not
a full chain, so *it* isn't cartoon-renderable), then renders the aligned
design's full protein as cartoon, its motif's `select_fixed_atoms` as
spheres, and the reference's ligand as licorice. When a motif residue's
named atoms don't survive in the folded design (LigandMPNN redesigned that
position's identity, since only RFDiffusion3's backbone generation — not
necessarily LigandMPNN's sequence design — treats it as fixed), it falls
back to that residue's backbone atoms so the position still reads as a hot
spot. See `examples/discontinuous_scaffolds_motif/` for a worked pipeline,
including the application-specific `find_best_fold.py` lookup step that
resolves a model name to its best (lowest `motif_rmsd`) completed run across
a campaign's `campaign_analysis.csv` files.

**`lib/peptide_align.py` is a fourth selection/alignment convention**,
alongside chain ids, `lib/ligand_select.py`'s atom-name ligand split, and
`lib/rfd3_motif_select.py`'s contig-based motif mapping: rigid-body Kabsch
superposition of one peptide chain's C-terminal residues onto another's,
matched purely by position (N -> C order) rather than residue identity or an
external contig map — `fit_cterminal_backbone` takes the last `n_residues`
of each chain (`cterminal_resi`) and Kabsch-fits their backbone (N, CA, C,
O) atoms (`backbone_coords`). This is the right tool when a short peptide
fragment (e.g. a PDZ domain's minimal conserved C-terminal binding motif,
used as a design template) and a longer peptide produced downstream share
that same C-terminal motif but otherwise differ in length or sequence, so no
sequence/structure alignment (`lib/rmsd.py`'s cealign) or explicit
correspondence map is needed. `aligned_pair_figure.py` (GENERATE) is built on
top of it: unlike `motif_superposition_figure.py`'s single overlaid panel, it
renders **two separate panels** from one PyMOL session — a `design` complex,
oriented with the standard `orient_long_axis_vertical` convention, and a
`reference` complex (e.g. the template structure a design was generated
from), whose coordinates are rigidly Kabsch-fit so its peptide's C-terminal
motif lands on `design`'s already-oriented peptide's equivalent residues.
Because orientation is baked into atom coordinates rather than being a
camera move (see below), and both panels are ray-traced from PyMOL's
untouched default view, the shared motif lands in the same place in both
renders even though the rest of each structure (domain sequence, peptide
length) differs freely — the two PNGs are meant to be tiled side by side
with `montage_figures.py` (a new `render_pair` `run_pipeline.py` step kind,
alongside `render`/`generate_each`, dispatches a script with two declared
outputs instead of one). The shared "load both objects, orient `design`,
Kabsch-fit `reference` into its frame" recipe lives in
`lib/peptide_align.py`'s `load_and_align_pair`, not duplicated per script.
`aligned_overlay_figure.py` is the same alignment rendered the other way —
**one** panel with both structures superposed (four independent colors:
reference domain/peptide, design domain/peptide, default cyan/green/red/purple)
rather than two tiled panels, for reading the domain and peptide change as a
direct overlay; it's a new `render_overlay` step kind (two inputs, one
output, unlike `render`'s one-input/one-output or `render_pair`'s
two-input/two-output). See `examples/pdz_design_vs_template/` for a worked
pipeline over the pdzbinder production campaign producing both views,
including the application-specific `resolve_template.py` lookup step that
maps a `find_structures_campaign.py` candidate's `group`/`campaign_dir` to
the *base* design's `prod_in/<base>_in/<id>.pdb` template (`campaign_base`,
a new `find_structures_campaign.py` helper that strips any `_subN`
iterative-refinement suffix — the specific stage that produced the winning
prediction is not, in general, the design lineage's actual starting point).

**`motif_superposition_figure.py` also supports a close-up, single-panel
mode** (`--ligand-representation surface` + `--orient-toward hotspot`),
added for `examples/discontinuous_scaffolds_motif_single/` without changing
any flag's default, so the five-panel pipeline above keeps rendering
identically. Two PyMOL gotchas surfaced building it, both handled in
`show_ligand()`:
- PyMOL silently auto-flags the small, disconnected hetero groups typical of
  an M-CSA-style minimal reference structure as `ignore` on load, which
  zeroes their surface area with no error — `cmd.flag("ignore", ..., "clear")`
  before `show surface` is required, or nothing renders.
- This close in, the camera is often inside the fold's own cartoon volume,
  so back-facing polygons are unavoidable; without `two_sided_lighting`,
  PyMOL renders those solid black instead of shaded-through.

`--orient-toward hotspot` does **not** reuse `orient_scene_vertically()`'s
"protein long axis vertical, twist toward a target" convention (that keeps
the *protein's* overall long axis vertical, which has nothing to do with
the motif's own local geometry and left the interface facing an arbitrary
direction). Instead it's `orient_look_down_on_plane()`: PCA over just the
hot-spot atoms gives a best-fit plane (smallest-variance eigenvector = the
plane's normal), which is rotated to face the camera (+Z) directly — i.e.
the camera looks straight down onto the interface — with the plane's own
largest-variance in-plane direction made vertical for a deterministic "up".
The normal's sign is picked so the protein's bulk (its CA atoms' centroid)
ends up behind the plane (-Z), not in front of it.

**Panning happens before any trim, always on the whole complex.** The
camera always pans/zooms to fit `reference or design` (everything actually
shown) with `--zoom-buffer` padding — never a tighter selection like just
the hot spot — so nothing rendered ends up cropped by the canvas edge (a
tight zoom straight to the motif could leave the wider cartoon extending
past the frame, which reads as the structure being cut off). Only *after*
that safe pan/zoom and the ray-trace does `build_figure()` reuse
`lib/imgtrim.py`'s `trim_to_content` (on by default, `--no-trim` to skip)
to crop the rendered PNG's excess background margin — a 2D crop of a
complete render, not a 3D camera move, so a close, content-filling
composition can never lose part of the structure.

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
quality. The tiling/scaling logic itself (`build_montage`, `scale_to_width`)
lives in `lib/montage.py`, shared with `assemble_panel_layout.py` rather than
each script carrying its own copy; both scripts are thin CLI wrappers over
it. Before tiling, `build_montage` trims each source image to its content
bounding box via `lib/imgtrim.py` (`trim_to_content`, `--no-trim` to skip) —
every GENERATE script rays-traces onto a fixed canvas with its own zoom
buffer, so raw PNGs otherwise carry inconsistent background margin that's
most visible when panels have different aspect ratios. `assemble_panel_layout.py`
handles a layout a uniform `rows x cols` grid can't express: a wide "problem"
panel spanning the full height of an adjacent `rows x cols` grid of "design"
panels (built by calling `lib/montage.py`'s `build_montage()` directly), the
left panel's width computed as a fraction of the *final* canvas.
