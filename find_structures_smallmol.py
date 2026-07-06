#!/usr/bin/env python3
"""
FIND step: discover completed small-molecule-binder designs produced by an
ImpressBasePipeline SmallMoleculeBindingPipeline run (see that project's
CLAUDE.md for the state machine this mirrors) — either a `prod` campaign
root (p1, p2, ... directly underneath) or a `nonadaptive_reference` one
(same p1, p2, ... layout, just without the adaptive-feedback loop).

Each pipeline campaign directory (`p<N>/`) contains a flat, numbered
sequence of HPC/local task directories: `<taskcount>_<taskname>/{in,out}/`.
A "completed design" is one that reached a fold-prediction (`alphafold`)
task with a computable mean pLDDT; this step also walks back one task
number (the fastrelax/filter_shape task number consistently immediately
precedes its alphafold task in every observed campaign here) to find the
FastRelax-packed, protein+ligand complex PDB — the AF2 output itself is
protein-only (no ligand: ColabFold doesn't co-fold small molecules), so the
FastRelax structure is the one actually worth rendering as "the design".

Usage:
    python3 find_structures_smallmol.py --campaign-root DIR --out candidates.json \
        [--campaign-label prod] [--chain-domain A] [--chain-peptide B]
"""
import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from manifest import write_manifest  # noqa: E402


def mean_plddt_and_source(alphafold_out_dir):
    """Best (highest mean pLDDT) rank among this task's *_scores_*.json files."""
    best_plddt, best_json = -1.0, None
    for sf in glob.glob(os.path.join(alphafold_out_dir, "*scores*.json")):
        try:
            with open(sf) as fh:
                arr = json.load(fh).get("plddt", [])
        except (OSError, json.JSONDecodeError):
            continue
        if not arr:
            continue
        mean = sum(arr) / len(arr)
        if mean > best_plddt:
            best_plddt = mean
            best_json = sf
    return (best_plddt, best_json) if best_json else (None, None)


def max_shape_complementarity(filter_shape_out_dir):
    sc_file = os.path.join(filter_shape_out_dir, "shape_complementarity_values.txt")
    if not os.path.exists(sc_file):
        return None
    max_sc = None
    with open(sc_file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    val = float(parts[1].split(": ")[-1])
                except ValueError:
                    continue
                max_sc = val if max_sc is None else max(max_sc, val)
    return max_sc


def fastrelax_complex_pdb(fastrelax_out_dir):
    pdbs = sorted(glob.glob(os.path.join(fastrelax_out_dir, "*.pdb")))
    return pdbs[0] if pdbs else None


def find_task_dirs(campaign_dir):
    """{task_num: {taskname: dirpath}} for every '<N>_<name>' subdir."""
    by_num = {}
    for entry in os.listdir(campaign_dir):
        full = os.path.join(campaign_dir, entry)
        if not os.path.isdir(full):
            continue
        m = re.match(r"^(\d+)_(.+)$", entry)
        if not m:
            continue
        by_num.setdefault(int(m.group(1)), {})[m.group(2)] = full
    return by_num


def discover_campaigns(campaign_root):
    return sorted(
        (e for e in os.listdir(campaign_root)
         if re.match(r"^p\d+$", e) and os.path.isdir(os.path.join(campaign_root, e))),
        key=lambda name: int(name[1:]),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--campaign-root", required=True, help="Directory containing p1, p2, ... campaign subdirectories")
    parser.add_argument("--campaign-label", default=None, help="Tag recorded on every candidate (e.g. 'prod', 'nonadaptive') for downstream grouping/labeling; default: basename of --campaign-root")
    parser.add_argument("--chain-domain", default="A", help="Chain id of the designed binder protein (default: A)")
    parser.add_argument("--chain-peptide", default="B", help="Chain id of the ligand hetero-group (default: B)")
    parser.add_argument("--out", required=True, help="Path to write the candidate manifest JSON")
    args = parser.parse_args()

    campaign_root = os.path.abspath(args.campaign_root)
    label = args.campaign_label or os.path.basename(campaign_root.rstrip("/"))
    campaigns = discover_campaigns(campaign_root)
    if not campaigns:
        sys.exit(f"No p<N> campaign directories found under {campaign_root}")

    candidates = []
    n_skipped = 0
    for campaign in campaigns:
        campaign_dir = os.path.join(campaign_root, campaign)
        by_num = find_task_dirs(campaign_dir)

        for num, types in sorted(by_num.items()):
            if "alphafold" not in types:
                continue
            plddt, af2_json = mean_plddt_and_source(os.path.join(types["alphafold"], "out"))
            if plddt is None:
                n_skipped += 1
                continue

            fr_num = num - 1
            fr_types = by_num.get(fr_num, {})
            if "fastrelax" not in fr_types:
                # fall back to the closest earlier fastrelax task in this campaign
                earlier = [n for n in by_num if n < num and "fastrelax" in by_num[n]]
                fr_num = max(earlier) if earlier else None
                fr_types = by_num.get(fr_num, {}) if fr_num is not None else {}

            if "fastrelax" not in fr_types:
                n_skipped += 1
                continue
            pdb_path = fastrelax_complex_pdb(os.path.join(fr_types["fastrelax"], "out"))
            if pdb_path is None or not os.path.exists(pdb_path):
                n_skipped += 1
                continue

            max_sc = max_shape_complementarity(os.path.join(fr_types["filter_shape"], "out")) \
                if "filter_shape" in fr_types else None

            candidates.append({
                "id": f"{label}_{campaign}_task{num}",
                "pdb_path": pdb_path,
                "campaign": campaign,
                "campaign_label": label,
                "task_num": num,
                "af2_plddt": plddt,
                "max_sc": max_sc,
                "chain_domain": args.chain_domain,
                "chain_peptide": args.chain_peptide,
            })

    write_manifest(args.out, candidates)
    print(
        f"Found {len(candidates)} completed design(s) across {len(campaigns)} campaign(s) under {campaign_root} "
        f"[{n_skipped} alphafold task(s) skipped: no plddt or no matching fastrelax pdb] -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
