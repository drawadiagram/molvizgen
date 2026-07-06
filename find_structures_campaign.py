#!/usr/bin/env python3
"""
FIND step: recursively discover every predicted-structure PDB file produced
by a production campaign laid out as

    <prod-root>/<group>/protein_binding/af_pipeline_outputs_multi/<campaign>/
        af/prediction/dimer_models/<id>/boltz_results_<id>/predictions/<id>/<id>_model_0.pdb

where <group> is a "p<N>-p<M>" batch directory and <campaign> is one of
"p<N>", "p<N>_sub1", "p<N>_sub1_sub2", "p<N>_sub1_sub2_sub3" (successive
optimization stages for the same target <id>).

For each structure file found, this step also looks for a co-located
"confidence_<id>_model_0.json" (Boltz confidence output) and, if present,
attaches its "confidence_score" to the candidate record — pure discovery,
no reduction: every stage variant of every target is emitted. A later FILTER
step (filter_best_score.py) is responsible for picking the best stage per
target.

This campaign's PDB files use multi-character, whitespace-tokenized chain
ids ("pdz"/"pep") rather than the single-character fixed-column chain ids
the PDB format actually reserves a column for. Naive fixed-column readers
(PyMOL included) truncate these to their first character, and "pdz"/"pep"
both truncate to "p" — silently colliding rather than erroring. This FIND
step normalizes such files (see lib/pdb_normalize.py) into cached copies
using standard single-character chain ids "A" (domain) / "B" (peptide), and
every candidate's manifest entry points at the normalized copy so downstream
FILTER/GENERATE steps can use plain "chain A"/"chain B" selectors exactly
like the current directory's data — no special-casing needed downstream.

Usage:
    python3 find_structures_campaign.py --prod-root DIR --out candidates.json \
        [--groups 'p1-p16,p17-p32'] \
        [--raw-chain-domain pdz] [--raw-chain-peptide pep] \
        [--normalized-dir DIR]
"""
import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from manifest import write_manifest  # noqa: E402
from pdb_normalize import detect_chain_tokens, normalize_pdb_chains  # noqa: E402

MODEL_GLOB = os.path.join(
    "{group}", "protein_binding", "af_pipeline_outputs_multi", "*",
    "af", "prediction", "dimer_models", "*", "**", "*_model_0.pdb",
)


def discover_groups(prod_root, groups_arg):
    if groups_arg:
        return [g.strip() for g in groups_arg.split(",") if g.strip()]
    entries = sorted(os.listdir(prod_root))
    return [e for e in entries if os.path.isdir(os.path.join(prod_root, e)) and re.match(r"^p\d+-p\d+$", e)]


def campaign_stage(campaign_dir_name, target_id):
    """Turn 'p15_sub1_sub2' (with target id '4e34') into the stage suffix
    'sub1_sub2', or '' for the bare 'p15' (no-sub) stage."""
    prefix_match = re.match(r"^p\d+", campaign_dir_name)
    if not prefix_match:
        return campaign_dir_name
    rest = campaign_dir_name[prefix_match.end():]
    return rest.lstrip("_")


def campaign_base(campaign_dir_name):
    """Turn 'p15_sub1_sub2' into the base design id 'p15' -- `campaign_stage`'s
    complement, stripping any iterative-refinement suffix rather than
    returning it. This is the id whose prod_in/<base>_in/ template seeded the
    *entire* design lineage (used e.g. by
    examples/pdz_design_vs_template/resolve_template.py to find the
    structure a design was originally generated from, as opposed to whatever
    input a later refinement round happened to start from)."""
    prefix_match = re.match(r"^p\d+", campaign_dir_name)
    return prefix_match.group() if prefix_match else campaign_dir_name


def load_confidence(pdb_path, target_id):
    conf_path = os.path.join(os.path.dirname(pdb_path), f"confidence_{target_id}_model_0.json")
    if not os.path.exists(conf_path):
        return None
    try:
        with open(conf_path) as f:
            return json.load(f).get("confidence_score")
    except (OSError, json.JSONDecodeError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prod-root", required=True, help="Root directory containing p*-p* campaign batches")
    parser.add_argument("--groups", default=None, help="Comma-separated batch dirs to include (default: all p*-p* under --prod-root)")
    parser.add_argument("--raw-chain-domain", default="pdz", help="Domain chain token as written in the raw campaign files (default: pdz)")
    parser.add_argument("--raw-chain-peptide", default="pep", help="Peptide chain token as written in the raw campaign files (default: pep)")
    parser.add_argument("--normalized-dir", default=None, help="Directory to write chain-normalized PDB copies (default: alongside --out, in a 'normalized_pdbs' subdir)")
    parser.add_argument("--out", required=True, help="Path to write the candidate manifest JSON")
    args = parser.parse_args()

    prod_root = os.path.abspath(args.prod_root)
    groups = discover_groups(prod_root, args.groups)
    if not groups:
        sys.exit(f"No p*-p* group directories found under {prod_root}")

    normalized_dir = args.normalized_dir or os.path.join(os.path.dirname(os.path.abspath(args.out)) or ".", "normalized_pdbs")
    chain_token_map = {args.raw_chain_domain: "A", args.raw_chain_peptide: "B"}

    candidates = []
    seen_paths = set()
    n_normalized = 0
    n_skipped = 0
    for group in groups:
        pattern = os.path.join(prod_root, MODEL_GLOB.format(group=group))
        for pdb_path in sorted(glob.glob(pattern, recursive=True)):
            if pdb_path in seen_paths:
                continue
            seen_paths.add(pdb_path)

            target_id = os.path.splitext(os.path.basename(pdb_path))[0]
            target_id = re.sub(r"_model_0$", "", target_id)

            m = re.search(r"af_pipeline_outputs_multi/([^/]+)/af/prediction/dimer_models/([^/]+)/", pdb_path)
            campaign_dir = m.group(1) if m else "unknown"
            model_target = m.group(2) if m else target_id
            stage = campaign_stage(campaign_dir, model_target)

            confidence_score = load_confidence(pdb_path, target_id)

            found_tokens = set(detect_chain_tokens(pdb_path))
            if found_tokens == {"A", "B"}:
                resolved_path = pdb_path
            elif found_tokens == set(chain_token_map):
                resolved_path = os.path.join(normalized_dir, f"{group}__{campaign_dir}__{target_id}.pdb")
                normalize_pdb_chains(pdb_path, resolved_path, chain_token_map)
                n_normalized += 1
            else:
                print(f"  SKIP (unexpected chain tokens {found_tokens}): {pdb_path}", file=sys.stderr)
                n_skipped += 1
                continue

            candidates.append({
                "id": target_id,
                "pdb_path": resolved_path,
                "source_pdb_path": pdb_path,
                "group": group,
                "campaign_dir": campaign_dir,
                "stage": stage,
                "confidence_score": confidence_score,
                "chain_domain": "A",
                "chain_peptide": "B",
            })

    write_manifest(args.out, candidates)
    n_targets = len({c["id"] for c in candidates})
    print(
        f"Found {len(candidates)} candidates ({n_targets} distinct targets) across {len(groups)} group(s) "
        f"[{n_normalized} chain-normalized, {n_skipped} skipped] -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
