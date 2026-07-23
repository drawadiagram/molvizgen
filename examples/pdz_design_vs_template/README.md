# pdz_design_vs_template

## Figure

Two views comparing the pdzbinder production campaign's single
highest-confidence design against the (usually much shorter-peptide)
template structure its lineage was originally generated from, both aligned
on the conserved C-terminal "EPEA" peptide motif:

- `panel_4k6y` — two separate, identically-oriented panels (domain cyan
  cartoon, peptide green cartoon in both) tiled side by side.
- `overlay_4k6y` — one combined panel with both structures superposed (four
  colors: template domain cyan, template peptide green, design domain red,
  design peptide yellow), plus solo exports of each structure alone since an
  overlay is easy to misread when the two cartoons occlude each other.

## Workflow

One application-specific lookup step, then two independent GENERATE (+
ASSEMBLE) branches over the same resolved pair:

1. **Resolve** (`resolve_template.sh`, prerequisite) — generic FIND
   (`find_structures_campaign.py`) + FILTER (`filter_top_n.py`) picks the
   single best-confidence production candidate; `resolve_template.py`
   (application-specific) looks up that candidate's originating `_in`
   template PDB by stripping any `_subN` iterative-refinement suffix from
   its `campaign_dir` (`find_structures_campaign.campaign_base`).
2. **GENERATE** (`render_pair` / `figures/aligned_pair_figure.py`) —
   `panel_4k6y`: two separate panels sharing a coordinate frame, the
   template's peptide Kabsch-fit onto the design's peptide's last 4
   residues (`lib/peptide_align.py`).
3. **ASSEMBLE** (`assemble`) — `row`: tile the two panels side by side.
4. **GENERATE** (`render_overlay` / `figures/aligned_overlay_figure.py`) —
   `overlay_4k6y`: the same alignment rendered as one superposed panel, plus
   solo exports of each structure.

## Run

```
./resolve_template.sh
cd examples/pdz_design_vs_template && python3 ../../run_pipeline.py design_vs_template_pipeline.yaml
```

or, both steps in sequence:

```
./run_design_vs_template.sh
```

## Files

| File | Purpose |
|---|---|
| `resolve_template.sh` | Prerequisite: FIND + FILTER + template lookup, writes `template_specs/<id>.json`. |
| `resolve_template.py` | Application-specific: resolves a winning candidate's `_in` template PDB. |
| `design_vs_template_pipeline.yaml` | The two-view pipeline (`render_pair` + `render_overlay`). |
| `run_design_vs_template.sh` | One-shot wrapper: `resolve_template.sh` then `run_pipeline.py`. |
