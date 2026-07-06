#!/usr/bin/env python3
"""
Application-specific FIND/lookup step for the discontinuous-scaffolds
motif-superposition figure (see ../../motif_superposition_figure.py and
../../lib/rfd3_motif_select.py). Not a general molvizgen FIND step: it knows
this one campaign's exact directory/CSV layout (an ImpressBasePipeline
DiscontinuousScaffoldsPipeline production run — see that project's
CLAUDE.md), the way its adaptive redesign branching numbers directories
(`disco_p<N>_0`, `..._R`, `..._R_R`, ...), and its `campaign_analysis.csv`
column schema (`model_name`, `motif_rmsd`, `run_dir`, `chai1_model_idx`).

For each requested model, this:
  1. Finds the `mcsa_41-N.json` RFD3 design-spec file (under
     --design-json-dir) that has a top-level entry for that model.
  2. Scans every `campaign_analysis.csv` under --prod-root (one per branch
     pipeline directory, root or redesign) for rows naming that model, and
     keeps the lowest-`motif_rmsd` row — mirroring the pipeline's own
     `check_fold_results()` best-fold selection (see discontinuous_scaffolds
     CLAUDE.md's "Three-stage execution flow" section), across however many
     redesign generations that model went through.
  3. Resolves that row's `run_dir`/`chai1_model_idx` to the actual Chai-1
     CIF path (`<branch_dir>/7_fold_pred/out/<run_dir>/prediction/pred.model_idx_<idx>.cif`),
     where `<branch_dir>` is the campaign_analysis.csv's own branch directory
     (three levels up: .../<branch_dir>/8_pipeline_analysis/out/campaign_analysis.csv).
  4. Writes one panel-spec JSON per model to --out-dir, consumed directly by
     motif_superposition_figure.py — every path in it is absolute.

Usage:
    python3 find_best_fold.py --prod-root /abs/path/to/prod \
        --design-json-dir /abs/path/to/scripts \
        --reference-pdb-dir /abs/path/to/scripts/mcsa_41 \
        --island-counts-csv /abs/path/to/scripts/island_counts.csv \
        --models M0110_1c0p M0157_1qh5 M0349_1e3v M0209_1lij M0078_1al6 \
        --out-dir /abs/path/to/out/panel_specs
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

DESIGN_JSON_RE = re.compile(r"^mcsa_41-\d+\.json$")


def find_design_json(design_json_dir, model_name):
    for path in sorted(glob.glob(os.path.join(design_json_dir, "mcsa_41-*.json"))):
        if not DESIGN_JSON_RE.match(os.path.basename(path)):
            continue
        with open(path) as f:
            spec = json.load(f)
        if model_name in spec:
            return os.path.abspath(path)
    return None


def find_best_row(prod_root, model_name):
    """Return (csv_path, row_dict) with the lowest motif_rmsd across every
    campaign_analysis.csv under prod_root that has a row for model_name."""
    best = None
    pattern = os.path.join(prod_root, "**", "campaign_analysis.csv")
    for csv_path in glob.glob(pattern, recursive=True):
        with open(csv_path, newline="") as f:
            for row in csv.DictReader(f):
                if row.get("model_name") != model_name:
                    continue
                try:
                    rmsd = float(row["motif_rmsd"])
                except (KeyError, ValueError):
                    continue
                if best is None or rmsd < best[1]:
                    best = (csv_path, rmsd, row)
    if best is None:
        return None, None
    csv_path, _rmsd, row = best
    return csv_path, row


def resolve_design_cif(csv_path, row):
    # csv_path = <branch_dir>/8_pipeline_analysis/out/campaign_analysis.csv
    branch_dir = os.path.dirname(os.path.dirname(os.path.dirname(csv_path)))
    return os.path.join(
        branch_dir, "7_fold_pred", "out", row["run_dir"], "prediction",
        f"pred.model_idx_{row['chai1_model_idx']}.cif",
    )


def load_island_counts(island_counts_csv):
    if not island_counts_csv:
        return {}
    with open(island_counts_csv, newline="") as f:
        return {row["INPUT_ID"]: int(row["RESIDUE_ISLAND_COUNT"]) for row in csv.DictReader(f) if row["INPUT_ID"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prod-root", required=True, help="Root directory to search recursively for campaign_analysis.csv files")
    parser.add_argument("--design-json-dir", required=True, help="Directory containing mcsa_41-N.json RFD3 design-spec files")
    parser.add_argument("--reference-pdb-dir", required=True, help="Directory containing reference PDBs, one per model (<model_name>.pdb)")
    parser.add_argument("--island-counts-csv", default=None, help="Optional island_counts.csv to attach RESIDUE_ISLAND_COUNT metadata")
    parser.add_argument("--rmsd-threshold", type=float, default=1.5, help="Recorded as pass/fail metadata only (default: 1.5)")
    parser.add_argument("--models", nargs="+", required=True, help="Model names to resolve (e.g. M0110_1c0p)")
    parser.add_argument("--out-dir", required=True, help="Directory to write one <model_name>.json panel spec per model")
    args = parser.parse_args()

    prod_root = os.path.abspath(args.prod_root)
    design_json_dir = os.path.abspath(args.design_json_dir)
    reference_pdb_dir = os.path.abspath(args.reference_pdb_dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    island_counts = load_island_counts(args.island_counts_csv)

    manifest = []
    for model_name in args.models:
        design_json = find_design_json(design_json_dir, model_name)
        if design_json is None:
            sys.exit(f"{model_name}: no mcsa_41-N.json under {design_json_dir} has an entry for this model")

        reference_pdb = os.path.join(reference_pdb_dir, f"{model_name}.pdb")
        if not os.path.exists(reference_pdb):
            sys.exit(f"{model_name}: reference PDB not found: {reference_pdb}")

        csv_path, row = find_best_row(prod_root, model_name)
        if row is None:
            sys.exit(f"{model_name}: no campaign_analysis.csv row found under {prod_root}")

        design_cif = resolve_design_cif(csv_path, row)
        if not os.path.exists(design_cif):
            sys.exit(f"{model_name}: resolved Chai-1 CIF not found: {design_cif} (from {csv_path})")

        motif_rmsd = float(row["motif_rmsd"])
        panel = {
            "model_name": model_name,
            "reference_pdb": reference_pdb,
            "design_json": design_json,
            "design_cif": design_cif,
            "motif_rmsd": motif_rmsd,
            "passes_rmsd_threshold": motif_rmsd < args.rmsd_threshold,
            "island_count": island_counts.get(model_name),
            "source_campaign_csv": csv_path,
            "run_dir": row["run_dir"],
            "chai1_model_idx": row["chai1_model_idx"],
        }
        out_path = os.path.join(out_dir, f"{model_name}.json")
        with open(out_path, "w") as f:
            json.dump(panel, f, indent=2)
        manifest.append(panel)
        print(f"{model_name}: motif_rmsd={motif_rmsd:.3f} island_count={panel['island_count']} -> {out_path}", file=sys.stderr)

    summary_path = os.path.join(out_dir, "_summary.json")
    with open(summary_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote {len(manifest)} panel spec(s) to {out_dir} (summary: {summary_path})", file=sys.stderr)


if __name__ == "__main__":
    main()
