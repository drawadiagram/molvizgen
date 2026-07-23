# diverse_figures

## Figure

A montage of 5 rendered PyMOL figures, one per pairwise-dissimilar PDZ-domain
structure selected (by chain-A Cα RMSD) from a directory of candidates —
the same diversity-screen idea as `../heatmap`, but with GENERATE + ASSEMBLE
stages added on top instead of stopping at a heatmap.

## Workflow

FIND → FILTER → GENERATE → ASSEMBLE, as a hand-written bash wrapper (not a
`run_pipeline.py` YAML) — useful as a second reference for how the
underlying scripts compose outside the YAML dispatcher, and for one-off
runs:

1. **FIND+FILTER** (`pdz_pairwise_rmsd.py`) — the older, directory-glob
   predecessor to `filter_diversity.py`: computes pairwise chain-A Cα RMSD
   over every PDB in a directory and greedily selects 5 dissimilar
   structures, no manifest involved.
2. **GENERATE** (`figures/pdz_figure.py`) — one auto-oriented PyMOL render
   per selected structure (fixed PDZ/peptide chain-A/chain-B convention,
   pink/lime-green coloring).
3. **ASSEMBLE** (`montage_figures.py`) — tile the 5 renders into one montage
   PNG.

## Run

```
./run_diverse_figures.sh [pdb-dir]
```

`pdb-dir` defaults to the same `new_decoys` directory `../heatmap` and
`../reference_experiment` use. No prerequisite steps.

## Files

| File | Purpose |
|---|---|
| `run_diverse_figures.sh` | The whole workflow: `pdz_pairwise_rmsd.py` → `figures/pdz_figure.py` (×5) → `montage_figures.py`. |
