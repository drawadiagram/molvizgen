# heatmap

## Figure

A pairwise Cα RMSD heatmap over a 5-structure, maximally-dissimilar subset
of a 70-structure candidate set — the simplest example in this repo, and a
good first one to read. No GENERATE/figure rendering at all: the output is
a matplotlib PNG, not a PyMOL render.

## Workflow

FIND → FILTER → PLOT, driven entirely by `rmsd_heatmap_pipeline.yaml`:

1. **FIND** (`find_flat`) — enumerate every `*.pdb` in a directory into a
   manifest.
2. **FILTER** (`filter_diversity`) — compute the full pairwise chain-domain
   Cα RMSD matrix (`lib/rmsd.py`, `cealign`) and greedily select the 5
   structures that are most mutually dissimilar.
3. **PLOT** (`plot_heatmap`) — render that RMSD matrix as an annotated
   heatmap PNG (`figures/plot_rmsd_heatmap.py`).

## Run

```
cd examples/heatmap && python3 ../../run_pipeline.py rmsd_heatmap_pipeline.yaml
```

No prerequisite steps — this is a self-contained FIND→FILTER→PLOT pipeline.

## Files

| File | Purpose |
|---|---|
| `rmsd_heatmap_pipeline.yaml` | The whole pipeline: `find_flat` → `filter_diversity` → `plot_heatmap`. |
