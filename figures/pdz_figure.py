#!/usr/bin/env python3
"""
Render a publication-style PyMOL figure of a single PDZ domain / peptide
complex: chain A (PDZ domain) as a pink cartoon, chain B (the PDZ-binding
motif peptide) as lime-green licorice, oriented so the domain's long axis
runs vertically, ray-traced with ambient occlusion and no shadows.

Usage:
    python3 pdz_figure.py <input.pdb> <output.png> \
        [--width 1800] [--height 2400] [--dpi 300] [--bg white]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from orient import orient_long_axis_vertical  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save  # noqa: E402

LIME_GREEN = [0x32 / 255, 0xCD / 255, 0x32 / 255]  # CSS "limegreen" (#32CD32)

PDZ_CHAIN = "chain A"
PEPTIDE_CHAIN = "chain B"


def build_figure(pdb_path, out_png, width, height, dpi, bg):
    obj = "complex"
    cmd.load(pdb_path, obj)
    cmd.remove(f"{obj} and hydro")

    orient_long_axis_vertical(
        f"{obj} and {PDZ_CHAIN} and name CA and polymer",
        f"{obj} and {PEPTIDE_CHAIN} and name CA and polymer",
        obj,
    )

    cmd.hide("everything", obj)

    # PDZ domain: cartoon, pink.
    cmd.show("cartoon", f"{obj} and {PDZ_CHAIN}")
    cmd.color("pink", f"{obj} and {PDZ_CHAIN}")

    # PBM peptide: bold licorice, lime green.
    cmd.show("sticks", f"{obj} and {PEPTIDE_CHAIN}")
    cmd.set_color("lime_green_custom", LIME_GREEN)
    cmd.set("stick_radius", 0.3, f"{obj} and {PEPTIDE_CHAIN}")
    cmd.set("stick_quality", 15)
    cmd.color("lime_green_custom", f"{obj} and {PEPTIDE_CHAIN}")

    apply_material_aoshiny()

    cmd.bg_color(bg)
    cmd.zoom(obj, buffer=5)

    ray_trace_and_save(out_png, width, height, dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb_path")
    parser.add_argument("output_png")
    add_render_flags(parser)
    args = parser.parse_args()

    build_figure(args.pdb_path, args.output_png, args.width, args.height, args.dpi, args.bg)
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
