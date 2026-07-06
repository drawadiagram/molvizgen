#!/usr/bin/env python3
"""
Application-specific lookup step for the pdz_design_vs_template example (see
../../aligned_pair_figure.py). Not a general molvizgen FIND step: FIND
(find_structures_campaign.py) and FILTER (filter_top_n.py/filter_best_score.py)
already resolve "which production complex(es)" generically -- this script
only knows the one extra, campaign-specific fact those general steps don't:
that every design lineage was itself generated starting from a template
structure at

    <prod-root>/<any-group>/protein_binding/prod_in/<base>_in/<id>.pdb

(a single PDB per base design, holding the pre-design domain sequence with
its own, often much shorter, copy of the bound peptide -- see the pdzbinder
project's own docs for why "prod_in" holds these; it's duplicated
byte-identically under every group directory, so which group's copy is read
doesn't matter).

A candidate's own `campaign_dir` (e.g. "p5_sub1_sub2" from
find_structures_campaign.py) names the *specific stage* that produced the
winning prediction, which may be several iterative LigandMPNN/RFDiffusion
refinement rounds ("_sub1", "_sub2", ...) downstream of the design's actual
origin -- those intermediate rounds' own "_in" input is just whatever the
previous round already produced (already carrying the full-length target
peptide), not the true originating template. The lineage's actual starting
point is always the *base* design id's template,
`prod_in/<base>_in/<id>.pdb`, where `<base>` is `campaign_dir` with any
"_subN" suffix stripped (find_structures_campaign.campaign_base, e.g.
"p5_sub1_sub2" -> "p5") -- usually a much shorter/minimal peptide fragment
(e.g. a PDZ target's bare "EPEA" C-terminal motif) than the full peptide the
production complex was ultimately generated against. Comparing the two is
the whole point of this example: visualizing how much the binder's peptide
interface grew over the course of the design lineage.

Given a manifest of already-selected production candidates (each already
carrying `group`/`campaign_dir`/`pdb_path` from find_structures_campaign.py),
this resolves each one's base-design `_in` template path and writes a small
`{id, group, campaign_dir, base_campaign, confidence_score, design_pdb,
reference_pdb}` JSON per candidate -- consumed directly by a
design_vs_template_pipeline.yaml `render_pair` step's `--design`/`--reference`
arguments, same convention as find_best_fold.py's panel-spec JSONs for
motif_superposition_figure.py.

Usage:
    python3 resolve_template.py --candidates best.json --prod-root /abs/path/to/prod \
        --out-dir /abs/path/to/out/template_specs
"""
import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from manifest import read_manifest  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from find_structures_campaign import campaign_base  # noqa: E402


def resolve_reference_pdb(prod_root, group, base):
    in_dir = os.path.join(prod_root, group, "protein_binding", "prod_in", f"{base}_in")
    matches = sorted(glob.glob(os.path.join(in_dir, "*.pdb")))
    if len(matches) != 1:
        sys.exit(f"expected exactly one *_in template PDB in {in_dir}, found {len(matches)}: {matches}")
    return matches[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidates", required=True, help="Manifest of already-selected production candidates (e.g. filter_top_n.py's output)")
    parser.add_argument("--prod-root", required=True, help="Root directory containing p*-p* campaign batches (same one passed to find_structures_campaign.py)")
    parser.add_argument("--out-dir", required=True, help="Directory to write one <id>.json template spec per candidate")
    args = parser.parse_args()

    prod_root = os.path.abspath(args.prod_root)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    candidates = read_manifest(args.candidates)
    if not candidates:
        sys.exit(f"{args.candidates} has no candidates")

    specs = []
    for c in candidates:
        base = campaign_base(c["campaign_dir"])
        reference_pdb = resolve_reference_pdb(prod_root, c["group"], base)
        spec = {
            "id": c["id"],
            "group": c["group"],
            "campaign_dir": c["campaign_dir"],
            "base_campaign": base,
            "confidence_score": c.get("confidence_score"),
            "design_pdb": c["pdb_path"],
            "reference_pdb": reference_pdb,
        }
        out_path = os.path.join(out_dir, f"{c['id']}.json")
        with open(out_path, "w") as f:
            json.dump(spec, f, indent=2)
        specs.append(spec)
        print(
            f"{c['id']} (group={c['group']}, campaign={c['campaign_dir']}, base={base}, confidence={c.get('confidence_score')}): "
            f"design={c['pdb_path']} reference={reference_pdb} -> {out_path}",
            file=sys.stderr,
        )

    summary_path = os.path.join(out_dir, "_summary.json")
    with open(summary_path, "w") as f:
        json.dump(specs, f, indent=2)
    print(f"Wrote {len(specs)} template spec(s) to {out_dir} (summary: {summary_path})", file=sys.stderr)


if __name__ == "__main__":
    main()
