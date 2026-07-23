# backbone_fold_overlay

## Figure

Up to seven images from one shared alignment: a de novo RFDiffusion3
backbone (one protomer) rendered alone, the downstream AlphaFold2/ColabFold
fold prediction of the *whole* oligomeric assembly it seeds rendered alone
(cartoon), two aligned overlap panels with representation/color swapped
between passes (`--alt-representation`, default surface), and three
independent fold-alone extras (surface, ribbon, and licorice-colored-by-
secondary-structure). The backbone/fold pair baked into the wrapper's
defaults (campaign `small_molecule_binding-2`, p2, task 73→78) was picked by
surveying every rfd3→alphafold lineage across p1-p8 for the best structural
agreement (RMSD 0.41 Å over 96/99 residues).

## Workflow

No FIND/FILTER step — this is a single hand-picked pair, not a filtered
set — just one GENERATE step producing all seven images:

1. **GENERATE** (`figures/backbone_fold_overlay_figure.py`) — the fold
   prediction is one long single-chain object (N copies of the protomer
   concatenated, continuous residue numbering) rather than N separate
   chains, so instead of a plain whole-structure `cealign` this uses
   `lib/oligomer_align.py` to slide a CA-count-sized window along the fold
   and align onto whichever window actually corresponds to the backbone
   (`find_best_repeat_window` + `align_onto_best_repeat`). Renders the
   backbone alone, the fold alone (cartoon), the two aligned overlap panels,
   and the fold-alone extras from that one alignment.

Bash-wrapper-only example (like `../diverse_figures`) — no `run_pipeline.py`
YAML, since a single one-off pair doesn't need a FIND/FILTER step around it.

## Run

```
./run_backbone_fold_overlay.sh [backbone.cif.gz] [fold.pdb] [out_dir]
```

All three arguments are optional and default to the picked
`small_molecule_binding-2` p2 task-73→78 pair described above.

## Files

| File | Purpose |
|---|---|
| `run_backbone_fold_overlay.sh` | The whole workflow: one call to `figures/backbone_fold_overlay_figure.py` producing all seven images. |
