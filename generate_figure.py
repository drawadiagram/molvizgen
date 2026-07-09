#!/usr/bin/env python3
"""
GENERATE step: render a publication-style PyMOL figure of a two-chain
domain/peptide complex — domain as a cartoon, peptide as bold licorice,
oriented so the domain's long axis runs vertically, ray-traced with ambient
occlusion and no shadows. Generalized from pdz_figure.py: the chain ids and
colors are CLI flags (with defaults matching pdz_figure.py's PDZ/peptide
look) instead of module constants, so the same rendering logic works for
applications that use different chain-id conventions (e.g. production
Boltz outputs using literal chain ids "pdz"/"pep" rather than "A"/"B").

Optionally, the peptide chain's coloring can instead follow the hot-spot/
rest split of an RFDiffusion3 binder-design spec (`--peptide-design-json`,
same spec/convention as `ligand_hotspot_figure.py`'s `lib/ligand_select.py`)
rather than a single flat `--peptide-color-hex` — this carries the
"basic problem" panel's red/yellow ligand coloring through into the bound
complex, so the same atoms read as targeted vs. not across every panel of a
figure.

Usage:
    python3 generate_figure.py <input.pdb> <output.png> \
        [--domain-chain A] [--peptide-chain B] \
        [--domain-color pink] [--peptide-color-hex 32CD32] \
        [--peptide-design-json design.json] [--peptide-design-task partial] \
        [--hotspot-color red] [--rest-color yellow] \
        [--width 1800] [--height 2400] [--dpi 300] [--bg white]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from ligand_select import atom_name_selection, load_binder_design_spec  # noqa: E402
from orient import orient_long_axis_vertical  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save  # noqa: E402


def hex_to_rgb01(hex_color):
    hex_color = hex_color.lstrip("#")
    return [int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def build_figure(pdb_path, out_png, domain_chain, peptide_chain, domain_color,
                  peptide_color_hex, width, height, dpi, bg,
                  peptide_design_json=None, peptide_design_task="partial",
                  hotspot_color="red", rest_color="yellow"):
    obj = "complex"
    cmd.load(pdb_path, obj)
    cmd.remove(f"{obj} and hydro")

    orient_long_axis_vertical(
        f"{obj} and chain {domain_chain} and name CA and polymer",
        f"{obj} and chain {peptide_chain} and name CA and polymer",
        obj,
    )

    cmd.hide("everything", obj)

    cmd.show("cartoon", f"{obj} and chain {domain_chain}")
    cmd.color(domain_color, f"{obj} and chain {domain_chain}")

    peptide_sel = f"{obj} and chain {peptide_chain}"
    cmd.show("sticks", peptide_sel)
    cmd.set("stick_radius", 0.3, peptide_sel)
    cmd.set("stick_quality", 15)

    if peptide_design_json:
        _, ligand, buried_atoms, exposed_atoms = load_binder_design_spec(
            peptide_design_json, peptide_design_task
        )
        ligand_sel = f"{peptide_sel} and resn {ligand}"
        hotspot_sel = atom_name_selection(ligand_sel, buried_atoms)
        rest_sel = atom_name_selection(ligand_sel, exposed_atoms) if exposed_atoms else \
            f"({ligand_sel}) and not ({hotspot_sel})"
        cmd.color(hotspot_color, hotspot_sel)
        cmd.color(rest_color, rest_sel)
    else:
        cmd.set_color("generate_figure_peptide_color", hex_to_rgb01(peptide_color_hex))
        cmd.color("generate_figure_peptide_color", peptide_sel)

    apply_material_aoshiny()

    cmd.bg_color(bg)
    cmd.zoom(obj, buffer=5)

    ray_trace_and_save(out_png, width, height, dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdb_path")
    parser.add_argument("output_png")
    parser.add_argument("--domain-chain", default="A", help="Chain id of the domain (default: A)")
    parser.add_argument("--peptide-chain", default="B", help="Chain id of the peptide (default: B)")
    parser.add_argument("--domain-color", default="pink", help="PyMOL color name for the domain cartoon (default: pink)")
    parser.add_argument("--peptide-color-hex", default="32CD32", help="Hex color (no #) for the peptide licorice (default: 32CD32 / limegreen); ignored when --peptide-design-json is given")
    parser.add_argument("--peptide-design-json", default=None,
                         help="Path to an RFDiffusion3 *_binder_design.json spec; when given, colors the peptide "
                              "chain's target-ligand atoms by hot-spot split (--hotspot-color/--rest-color) instead "
                              "of a single flat --peptide-color-hex")
    parser.add_argument("--peptide-design-task", default="partial", choices=["buried", "partial"],
                         help="Which task's atom split to use from --peptide-design-json (default: partial)")
    parser.add_argument("--hotspot-color", default="red", help="PyMOL color for the targeted/buried peptide atoms when --peptide-design-json is given (default: red)")
    parser.add_argument("--rest-color", default="yellow", help="PyMOL color for the rest of the peptide when --peptide-design-json is given (default: yellow)")
    add_render_flags(parser)
    args = parser.parse_args()

    build_figure(
        args.pdb_path, args.output_png, args.domain_chain, args.peptide_chain,
        args.domain_color, args.peptide_color_hex, args.width, args.height, args.dpi, args.bg,
        peptide_design_json=args.peptide_design_json, peptide_design_task=args.peptide_design_task,
        hotspot_color=args.hotspot_color, rest_color=args.rest_color,
    )
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
