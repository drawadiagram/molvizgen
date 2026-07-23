#!/usr/bin/env python3
"""
Application-specific FIND/lookup step for a two-panel "design progression"
figure (see ../../figures/motif_superposition_figure.py and
../../lib/rfd3_motif_select.py's `anchor_chai_positions`). Not a general
molvizgen FIND kind: it knows this one campaign's exact directory layout (an
ImpressBasePipeline DiscontinuousScaffoldsPipeline production run — see that
project's CLAUDE.md), its adaptive fold-redesign branch naming
(`disco_p<N>_0` for the root pipeline, `..._R`, `..._R_R`, ... for
successive redesign generations of a model that failed `rmsd_threshold`),
and its `campaign_analysis.csv` column schema (`model_name`, `motif_rmsd`,
`run_dir`, `chai1_model_idx`, `anchor_residues`).

A single input motif is typically processed by *several independent*
pipeline runs (one per top-level campaign batch directory under
--prod-root, e.g. `p9-p16/` and `p11-p31r/` both contain a `disco_p14_*`
redesign lineage for model M0157_1qh5, from two separate campaigns) — not
every run's redesign lineage actually reaches a passing model before
hitting the pipeline's MAX_REDESIGN_DEPTH/MIN_RMSD_IMPROVEMENT stopping
conditions. This script:

  1. Scans every `<prod_root>/*/discontinuous_scaffolds` batch directory's
     `scripts/mcsa_41-N.json` files for one with a top-level entry for
     --model, to learn that batch's protein index N and root design-spec
     path.
  2. Within each such batch, walks the `disco_p<N>_0`, `disco_p<N>_R`,
     `disco_p<N>_R_R`, ... redesign lineage (in order), reading each
     generation's `8_pipeline_analysis/out/campaign_analysis.csv` (skipping
     any generation whose analysis hasn't completed) and keeping the
     lowest-`motif_rmsd` row per generation — mirroring the pipeline's own
     `check_fold_results()` best-fold selection (see discontinuous_scaffolds
     CLAUDE.md's "Three-stage execution flow" section).
  3. Picks the one batch/lineage whose deepest completed generation passes
     --rmsd-threshold and required at least one redesign (a lineage whose
     very first/root generation already passes has no "anchoring" to show
     in a second panel) — the pipeline run "that ends with a passing model".
     Use --batch to force a specific one instead of auto-selecting.
  4. Computes the *anchor* chai-sequence positions from the root
     generation's own `anchor_residues` field (the residues identified as
     well-predicted enough to hold fixed in a subsequent redesign) mapped
     through the root's own contig (`lib/rfd3_motif_select.anchor_chai_positions`).
     Because that position arithmetic is invariant across redesign
     generations, the *same* position set is reused to highlight "the
     portion that was subject to anchoring" in the final generation's own
     motif too.
  5. Writes two panel-spec JSONs — `<model>_initial.json` (the root
     generation) and `<model>_final.json` (the passing generation) — each
     consumed directly by motif_superposition_figure.py. Every path in them
     is absolute.

Usage:
    python3 find_anchor_progression.py --prod-root /abs/path/to/prod \
        --model M0157_1qh5 --rmsd-threshold 1.5 \
        --out-dir /abs/path/to/out/panel_specs
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "lib"))
from rfd3_motif_select import load_design_spec, anchor_chai_positions  # noqa: E402

DESIGN_JSON_RE = re.compile(r"^mcsa_41-(\d+)\.json$")
GENERATION_RE = re.compile(r"^disco_p(\d+)((?:_R)+|_0)$")


def batch_label(batch_dir):
    """Human-readable batch name for a '<prod_root>/<batch>/discontinuous_scaffolds' path, e.g. 'p11-p31r'."""
    return os.path.basename(os.path.dirname(batch_dir))


def discover_batches(prod_root):
    return sorted(
        d for d in glob.glob(os.path.join(prod_root, "*", "discontinuous_scaffolds"))
        if os.path.isdir(d)
    )


def find_root_design(batch_dir, model_name):
    """Return (protein_idx, design_json_path) for the mcsa_41-N.json under
    batch_dir/scripts that has a top-level entry for model_name, or
    (None, None) if this batch never processed that model."""
    scripts_dir = os.path.join(batch_dir, "scripts")
    for path in sorted(glob.glob(os.path.join(scripts_dir, "mcsa_41-*.json"))):
        m = DESIGN_JSON_RE.match(os.path.basename(path))
        if not m:
            continue
        with open(path) as f:
            spec = json.load(f)
        if model_name in spec:
            return int(m.group(1)), os.path.abspath(path)
    return None, None


def lineage_generations(batch_dir, protein_idx):
    """Return [(depth, branch_name), ...] sorted root-first for every
    disco_p<protein_idx>_(0|R+) directory directly under batch_dir (ignoring
    unrelated backbone/sequence adaptive branches like disco_p<N>_1)."""
    generations = []
    for entry in os.listdir(batch_dir):
        m = GENERATION_RE.match(entry)
        if not m or int(m.group(1)) != protein_idx:
            continue
        suffix = m.group(2)
        depth = 0 if suffix == "_0" else suffix.count("R")
        generations.append((depth, entry))
    return sorted(generations)


def best_row(csv_path, model_name):
    """Lowest-motif_rmsd row for model_name in one campaign_analysis.csv, or
    None if no row matches."""
    best = None
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("model_name") != model_name:
                continue
            try:
                rmsd = float(row["motif_rmsd"])
            except (KeyError, ValueError):
                continue
            if best is None or rmsd < best[0]:
                best = (rmsd, row)
    return best[1] if best else None


def resolve_design_cif(branch_dir, row):
    return os.path.join(
        branch_dir, "7_fold_pred", "out", row["run_dir"], "prediction",
        f"pred.model_idx_{row['chai1_model_idx']}.cif",
    )


def evaluate_batch(batch_dir, model_name, rmsd_threshold):
    """Return a dict describing this batch's disco_p<N> redesign lineage for
    model_name (generations, final generation, pass/fail), or None if this
    batch never produced completed analysis for the model at all."""
    protein_idx, root_design_json = find_root_design(batch_dir, model_name)
    if protein_idx is None:
        return None

    generations = []
    for depth, branch_name in lineage_generations(batch_dir, protein_idx):
        branch_dir = os.path.join(batch_dir, branch_name)
        csv_path = os.path.join(branch_dir, "8_pipeline_analysis", "out", "campaign_analysis.csv")
        if not os.path.exists(csv_path):
            continue  # generation still running / never reached the fold-analysis stage
        row = best_row(csv_path, model_name)
        if row is None:
            continue
        generations.append({
            "depth": depth,
            "branch_name": branch_name,
            "branch_dir": branch_dir,
            "csv_path": csv_path,
            "row": row,
            "motif_rmsd": float(row["motif_rmsd"]),
        })
    if not generations:
        return None

    final = generations[-1]
    return {
        "batch_dir": batch_dir,
        "protein_idx": protein_idx,
        "root_design_json": root_design_json,
        "generations": generations,
        "final": final,
        "passes": final["motif_rmsd"] < rmsd_threshold,
        "redesigned": final["depth"] > 0,
    }


def choose_lineage(candidates, forced_batch):
    if forced_batch is not None:
        for c in candidates:
            if batch_label(c["batch_dir"]) == forced_batch:
                return c
        sys.exit(f"--batch {forced_batch!r} not among batches with data: "
                  f"{[batch_label(c['batch_dir']) for c in candidates]}")

    passing_redesigned = [c for c in candidates if c["passes"] and c["redesigned"]]
    if len(passing_redesigned) == 1:
        return passing_redesigned[0]
    if len(passing_redesigned) > 1:
        names = [batch_label(c["batch_dir"]) for c in passing_redesigned]
        sys.exit(f"Multiple pipeline runs end in a passing, redesigned model: {names} "
                  f"— pick one with --batch")

    passing_any = [c for c in candidates if c["passes"]]
    if len(passing_any) == 1:
        return passing_any[0]

    lines = [
        f"  {batch_label(c['batch_dir'])}: "
        f"final={c['final']['branch_name']} motif_rmsd={c['final']['motif_rmsd']:.3f} "
        f"passes={c['passes']}"
        for c in candidates
    ]
    sys.exit("No single pipeline run unambiguously ends in a passing model:\n" + "\n".join(lines))


def build_panel(model_name, generation, design_json, reference_pdb, anchor_positions,
                 generation_label, rmsd_threshold, batch_dir):
    design_cif = resolve_design_cif(generation["branch_dir"], generation["row"])
    if not os.path.exists(design_cif):
        sys.exit(f"{model_name}: resolved Chai-1 CIF not found: {design_cif}")
    return {
        "model_name": model_name,
        "generation_label": generation_label,
        "generation_branch": generation["branch_name"],
        "generation_depth": generation["depth"],
        "reference_pdb": reference_pdb,
        "design_json": design_json,
        "design_cif": design_cif,
        "motif_rmsd": generation["motif_rmsd"],
        "passes_rmsd_threshold": generation["motif_rmsd"] < rmsd_threshold,
        "anchor_positions": sorted(anchor_positions),
        "anchor_residues_source": generation["row"].get("anchor_residues", ""),
        "source_campaign_csv": generation["csv_path"],
        "run_dir": generation["row"]["run_dir"],
        "chai1_model_idx": generation["row"]["chai1_model_idx"],
        "batch": batch_label(batch_dir),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prod-root", required=True, help="discontinuous_scaffolds prod root (contains p1-p8/, p9-p16/, p11-p31r/, ...)")
    parser.add_argument("--model", required=True, help="Model name to resolve (e.g. M0157_1qh5)")
    parser.add_argument("--rmsd-threshold", type=float, default=1.5, help="Passing threshold for motif_rmsd (default: 1.5)")
    parser.add_argument("--batch", default=None, help="Force a specific batch directory name (e.g. p11-p31r) instead of auto-selecting the one whose lineage ends in a passing, redesigned model")
    parser.add_argument("--out-dir", required=True, help="Directory to write <model>_initial.json / <model>_final.json panel specs")
    args = parser.parse_args()

    prod_root = os.path.abspath(args.prod_root)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    candidates = []
    for batch_dir in discover_batches(prod_root):
        result = evaluate_batch(batch_dir, args.model, args.rmsd_threshold)
        if result is not None:
            candidates.append(result)
    if not candidates:
        sys.exit(f"No campaign_analysis.csv rows found for {args.model} under {prod_root}")

    chosen = choose_lineage(candidates, args.batch)
    root = chosen["generations"][0]
    final = chosen["final"]

    scripts_dir = os.path.join(chosen["batch_dir"], "scripts")
    reference_pdb = os.path.join(scripts_dir, "mcsa_41", f"{args.model}.pdb")
    if not os.path.exists(reference_pdb):
        sys.exit(f"{args.model}: reference PDB not found: {reference_pdb}")

    root_entry = load_design_spec(chosen["root_design_json"], args.model)
    anchor_positions = anchor_chai_positions(root["row"].get("anchor_residues", ""), root_entry["contig"])

    # Every generation's design_json here is the *root* mcsa_41-N.json, not
    # a redesign generation's own redesign.json — a redesign's contig/
    # select_fixed_atoms keys are renumbered relative to the true reference
    # PDB (see create_redesign.py), so they can't be used as
    # motif_superposition_figure.py's `reference and chain ... resi ...`
    # selector. The root contig's chai-sequence-position arithmetic is
    # invariant across redesign generations (identical free-residue run
    # lengths, only the fixed-residue *labels* change — see this project's
    # CLAUDE.md and lib/rfd3_motif_select.py's module docstring), so reusing
    # it against a later generation's folded design_cif still resolves each
    # true motif residue to its correct position in that fold.
    initial_panel = build_panel(
        args.model, root, chosen["root_design_json"], reference_pdb, anchor_positions,
        "initial", args.rmsd_threshold, chosen["batch_dir"],
    )
    final_panel = build_panel(
        args.model, final, chosen["root_design_json"], reference_pdb, anchor_positions,
        "final", args.rmsd_threshold, chosen["batch_dir"],
    )

    initial_path = os.path.join(out_dir, f"{args.model}_initial.json")
    final_path = os.path.join(out_dir, f"{args.model}_final.json")
    with open(initial_path, "w") as f:
        json.dump(initial_panel, f, indent=2)
    with open(final_path, "w") as f:
        json.dump(final_panel, f, indent=2)

    print(
        f"{args.model}: batch={batch_label(chosen['batch_dir'])} "
        f"initial={root['branch_name']} (motif_rmsd={root['motif_rmsd']:.3f}) -> {initial_path}\n"
        f"{args.model}: batch={batch_label(chosen['batch_dir'])} "
        f"final={final['branch_name']} (motif_rmsd={final['motif_rmsd']:.3f}) -> {final_path}\n"
        f"{args.model}: anchor_positions (chai sequence positions) = {sorted(anchor_positions)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
