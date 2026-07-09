#!/usr/bin/env python3
"""
Generic runner for a YAML-described SELECT|GENERATE|ASSEMBLE pipeline.

A pipeline is a list of named steps:

    steps:
      - name: ref_select
        kind: filter_diversity
        args: {in: ${ref_find.manifest}, chain_field: chain_domain, n_select: 5}

Each step's `kind` dispatches to one of the existing CLI scripts in this
directory (find_structures_flat.py, find_structures_campaign.py,
filter_best_score.py, filter_diversity.py, plot_rmsd_heatmap.py,
montage_figures.py) or to the generate_each pseudo-step, which loops a
GENERATE script over a manifest.
`args` keys map 1:1 onto the underlying script's CLI flags
(`n_select` -> `--n-select`); a later step can reference any earlier step's
declared outputs with `${step_name.field}` (e.g. `${ref_find.manifest}`,
`${ref_figures.pngs}`, `${ref_montage.image}`) — this form must be the
*entire* value of a key and resolves to whatever type that field holds.

Data vs. transformation: a pipeline YAML declares only the transformation
(steps, filters, colors, layout) plus two small pieces of bookkeeping —
`figure_name` (a required descriptive label, used as the output subfolder
name) and an optional top-level `data:` block of named input-directory
defaults. Actual data locations are an invocation-time concern:

    figure_name: my_comparison
    out_root: /path/to/output/base       # optional; default cwd "."
    data:
      prod: /path/to/prod/experiment-1   # optional defaults

    steps:
      - name: prod_find
        args: {campaign_root: "${data.prod}"}

`${data.NAME}` resolves against the merged `data:` block + `--root`
overrides, and — unlike `${step.field}` — may be embedded inside a larger
string (e.g. `"${data.prod}/subdir/file.json"`), since data roots are
always plain strings.

Usage:
    python3 run_pipeline.py pipeline.yaml
    python3 run_pipeline.py pipeline.yaml --root prod=/path/to/experiment-2
    python3 run_pipeline.py pipeline.yaml --root name=/path --out-root /other/base
"""
import argparse
import os
import re
import subprocess
import sys

import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REF_PATTERN = re.compile(r"^\$\{([\w.]+)\}$")
DATA_REF_PATTERN = re.compile(r"\$\{data\.([\w-]+)\}")

sys.path.insert(0, os.path.join(SCRIPT_DIR, "lib"))
from manifest import read_manifest  # noqa: E402


def script_path(name):
    return os.path.join(SCRIPT_DIR, name)


def run(cmd):
    print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def flag(key):
    return "--" + key.replace("_", "-")


def flags_from_args(args, skip=()):
    """Generic dict-of-args -> CLI flags mapper: {'n_select': 5} -> ['--n-select', '5'].
    Booleans become bare store_true-style flags (True -> included, False -> omitted).
    Skips any key named in `skip` (handled specially by the caller)."""
    cmd = []
    for key, value in args.items():
        if key in skip or value is None:
            continue
        if isinstance(value, bool):
            if value:
                cmd.append(flag(key))
        else:
            cmd += [flag(key), str(value)]
    return cmd


def resolve_data(name, roots, context_value):
    if name not in roots:
        sys.exit(
            f"Reference to undeclared data root {name!r} in '{context_value}' "
            f"(declared: {sorted(roots) or 'none'}; add it to the YAML's 'data:' "
            f"block or pass --root {name}=PATH)"
        )
    return roots[name]


def resolve(value, ctx, roots):
    """Recursively substitute ${step_name.field} and ${data.name} references.
    `ctx` is a dict of {step_name: {field: value}} for prior steps' declared
    outputs; `roots` is a dict of {name: path} for declared/supplied data
    roots. A value that is *entirely* one ${...} placeholder resolves to
    whatever type that field holds (str or list); placeholders embedded in a
    list are resolved element-wise. ${data.name} may additionally appear
    embedded inside a larger string (e.g. "${data.prod}/subdir/file.json"),
    since data roots are always plain strings — ${step_name.field} may not."""
    if isinstance(value, str):
        m = REF_PATTERN.match(value)
        if m:
            step_name, _, field = m.group(1).partition(".")
            if step_name == "data":
                return resolve_data(field, roots, value)
            if step_name not in ctx:
                sys.exit(f"Reference to unknown step {step_name!r} in '{value}'")
            if field not in ctx[step_name]:
                sys.exit(f"Step {step_name!r} has no output {field!r} (has: {list(ctx[step_name])})")
            return ctx[step_name][field]
        if DATA_REF_PATTERN.search(value):
            return DATA_REF_PATTERN.sub(lambda mm: resolve_data(mm.group(1), roots, value), value)
        return value
    if isinstance(value, list):
        return [resolve(v, ctx, roots) for v in value]
    if isinstance(value, dict):
        return {k: resolve(v, ctx, roots) for k, v in value.items()}
    return value


def step_out_dir(base_out_dir, name):
    d = os.path.join(base_out_dir, name)
    os.makedirs(d, exist_ok=True)
    return d


# --- step kind handlers -----------------------------------------------------
# Each handler receives (name, args, base_out_dir) with args already
# ${...}-resolved, runs the underlying tool(s), and returns a dict of
# declared outputs available to later steps as ${name.<key>}.

def handle_find_flat(name, args, base_out_dir):
    manifest = os.path.join(step_out_dir(base_out_dir, name), "candidates.json")
    cmd = [sys.executable, script_path("find_structures_flat.py")]
    cmd += flags_from_args(args)
    cmd += ["--out", manifest]
    run(cmd)
    return {"manifest": manifest}


def handle_find_campaign(name, args, base_out_dir):
    out_dir = step_out_dir(base_out_dir, name)
    manifest = os.path.join(out_dir, "candidates.json")
    cmd = [sys.executable, script_path("find_structures_campaign.py")]
    cmd += flags_from_args(args, skip=("normalized_dir",))
    cmd += ["--normalized-dir", args.get("normalized_dir", os.path.join(out_dir, "normalized_pdbs"))]
    cmd += ["--out", manifest]
    run(cmd)
    return {"manifest": manifest}


def handle_find_smallmol(name, args, base_out_dir):
    manifest = os.path.join(step_out_dir(base_out_dir, name), "candidates.json")
    cmd = [sys.executable, script_path("find_structures_smallmol.py")]
    cmd += flags_from_args(args)
    cmd += ["--out", manifest]
    run(cmd)
    return {"manifest": manifest}


def handle_filter_top_n(name, args, base_out_dir):
    manifest = os.path.join(step_out_dir(base_out_dir, name), "candidates.json")
    cmd = [sys.executable, script_path("filter_top_n.py")]
    cmd += flags_from_args(args, skip=("in",))
    cmd += ["--in", args["in"]]
    cmd += ["--out", manifest]
    run(cmd)
    return {"manifest": manifest}


def handle_filter_best_score(name, args, base_out_dir):
    manifest = os.path.join(step_out_dir(base_out_dir, name), "candidates.json")
    cmd = [sys.executable, script_path("filter_best_score.py")]
    cmd += flags_from_args(args, skip=("in",))
    cmd += ["--in", args["in"]]
    cmd += ["--out", manifest]
    run(cmd)
    return {"manifest": manifest}


def handle_filter_diversity(name, args, base_out_dir):
    out_dir = step_out_dir(base_out_dir, name)
    cmd = [sys.executable, script_path("filter_diversity.py")]
    cmd += flags_from_args(args, skip=("in",))
    if "in" in args:
        cmd += ["--in", args["in"]]
    cmd += ["--out-dir", out_dir]
    run(cmd)
    manifest = os.path.join(out_dir, "selection.json")
    matrix = os.path.join(out_dir, "rmsd_matrix.csv")
    return {"manifest": manifest, "out_dir": out_dir, "matrix": matrix}


def handle_plot_heatmap(name, args, base_out_dir):
    out_dir = step_out_dir(base_out_dir, name)
    image = os.path.join(out_dir, args.get("out", "heatmap.png"))
    cmd = [sys.executable, script_path("plot_rmsd_heatmap.py")]
    cmd += flags_from_args(args, skip=("out",))
    cmd += ["--out", image]
    run(cmd)
    return {"image": image}


def handle_generate_each(name, args, base_out_dir):
    out_dir = step_out_dir(base_out_dir, name)
    gen_script = script_path(args.get("script", "pdz_figure.py"))
    extra = flags_from_args(args, skip=("selection", "script", "out_dir"))

    candidates = read_manifest(args["selection"])
    pngs = []
    for c in candidates:
        out_png = os.path.join(out_dir, f"{c['id']}_complex.png")
        cmd = [sys.executable, gen_script, c["pdb_path"], out_png] + extra
        run(cmd)
        pngs.append(out_png)
    return {"pngs": pngs, "out_dir": out_dir}


def handle_assemble(name, args, base_out_dir):
    out_dir = step_out_dir(base_out_dir, name)
    image = os.path.join(out_dir, args.get("out", "montage.png"))
    images = args["images"]
    if isinstance(images, str):
        images = [images]
    cmd = [sys.executable, script_path("montage_figures.py")]
    cmd += images
    cmd += flags_from_args(args, skip=("images", "out"))
    cmd += ["--out", image]
    run(cmd)
    return {"image": image}


def handle_render(name, args, base_out_dir):
    """One-off GENERATE call for a script that takes a single (input, output)
    pair rather than looping over a manifest (e.g. ligand_hotspot_figure.py,
    which renders one design-spec JSON, not a set of candidate structures)."""
    out_dir = step_out_dir(base_out_dir, name)
    image = os.path.join(out_dir, args.get("out", f"{name}.png"))
    cmd = [sys.executable, script_path(args["script"])]
    cmd += [args["input"], image]
    cmd += flags_from_args(args, skip=("script", "input", "out"))
    run(cmd)
    return {"image": image}


def handle_render_pair(name, args, base_out_dir):
    """GENERATE call for a script that renders two coupled panels from one
    shared alignment in a single PyMOL session (aligned_pair_figure.py)
    rather than a single (input, output) pair (handle_render) or a
    per-candidate loop (handle_generate_each)."""
    out_dir = step_out_dir(base_out_dir, name)
    reference_image = os.path.join(out_dir, args.get("out_reference", "reference.png"))
    design_image = os.path.join(out_dir, args.get("out_design", "design.png"))
    cmd = [sys.executable, script_path(args["script"])]
    cmd += [args["reference"], args["design"], reference_image, design_image]
    cmd += flags_from_args(args, skip=("script", "reference", "design", "out_reference", "out_design"))
    run(cmd)
    return {"reference_image": reference_image, "design_image": design_image}


def handle_render_overlay(name, args, base_out_dir):
    """GENERATE call for a script that overlays two structures into a single
    combined panel from a shared alignment (aligned_overlay_figure.py) —
    two positional inputs (reference, design) and one primary output, unlike
    handle_render (one input, one output) or handle_render_pair (two
    inputs, two outputs). Also resolves the optional --out-reference/
    --out-design solo-view exports every overlay-style script should offer
    (see aligned_overlay_figure.py) into this step's out_dir, declaring them
    as reference_image/design_image when requested."""
    out_dir = step_out_dir(base_out_dir, name)
    image = os.path.join(out_dir, args.get("out", f"{name}.png"))
    skip = ["script", "reference", "design", "out", "out_reference", "out_design"]
    cmd = [sys.executable, script_path(args["script"])]
    cmd += [args["reference"], args["design"], image]

    outputs = {"image": image}
    for key, field in (("out_reference", "reference_image"), ("out_design", "design_image")):
        if key in args:
            path = os.path.join(out_dir, args[key])
            cmd += [flag(key), path]
            outputs[field] = path

    cmd += flags_from_args(args, skip=skip)
    run(cmd)
    return outputs


def handle_assemble_panel_layout(name, args, base_out_dir):
    out_dir = step_out_dir(base_out_dir, name)
    image = os.path.join(out_dir, args.get("out", "panel_layout.png"))
    cmd = [sys.executable, script_path("assemble_panel_layout.py")]
    cmd += flags_from_args(args, skip=("left", "right", "out"))
    cmd += ["--left", args["left"]]
    cmd += ["--right"] + args["right"]
    cmd += ["--out", image]
    run(cmd)
    return {"image": image}


HANDLERS = {
    "find_flat": handle_find_flat,
    "find_campaign": handle_find_campaign,
    "find_smallmol": handle_find_smallmol,
    "filter_best_score": handle_filter_best_score,
    "filter_top_n": handle_filter_top_n,
    "filter_diversity": handle_filter_diversity,
    "plot_heatmap": handle_plot_heatmap,
    "generate_each": handle_generate_each,
    "render": handle_render,
    "render_pair": handle_render_pair,
    "render_overlay": handle_render_overlay,
    "assemble": handle_assemble,
    "assemble_panel_layout": handle_assemble_panel_layout,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="Path to a pipeline YAML file")
    parser.add_argument(
        "--root", action="append", default=[], metavar="[NAME=]PATH",
        help="Supply/override a named data root referenced as ${data.NAME}. "
             "Bare PATH (no 'NAME=') sets the root named 'default'. Repeatable.",
    )
    parser.add_argument(
        "--out-root", default=None,
        help="Override the YAML's out_root default. Output lands at <out-root>/<figure_name>/.",
    )
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    figure_name = config.get("figure_name")
    if not figure_name:
        sys.exit("Pipeline YAML must declare a top-level 'figure_name'")

    if any(step.get("name") == "data" for step in config["steps"]):
        sys.exit("Step name 'data' is reserved for the ${data.NAME} namespace; rename this step.")

    roots = {
        name: os.path.abspath(path)
        for name, path in (config.get("data") or {}).items()
        if path is not None
    }
    for raw in args.root:
        name, sep, path = raw.partition("=")
        if not sep:
            name, path = "default", raw
        roots[name] = os.path.abspath(path)

    out_root = args.out_root or config.get("out_root", ".")
    base_out_dir = os.path.abspath(os.path.join(out_root, figure_name))
    os.makedirs(base_out_dir, exist_ok=True)

    ctx = {}
    for step in config["steps"]:
        name = step["name"]
        kind = step["kind"]
        if kind not in HANDLERS:
            sys.exit(f"Unknown step kind {kind!r} in step {name!r} (known: {sorted(HANDLERS)})")

        print(f"== Step '{name}' ({kind}) ==", file=sys.stderr)
        step_args = resolve(step.get("args", {}), ctx, roots)
        ctx[name] = HANDLERS[kind](name, step_args, base_out_dir)

    print("", file=sys.stderr)
    print("Done. Step outputs:", file=sys.stderr)
    for name, outputs in ctx.items():
        for field, value in outputs.items():
            print(f"  {name}.{field} = {value}", file=sys.stderr)


if __name__ == "__main__":
    main()
