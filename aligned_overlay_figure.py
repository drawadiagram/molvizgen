#!/usr/bin/env python3
"""
GENERATE step: render a *single* overlaid panel comparing a design-template
structure ("reference") and the resulting production complex ("design"),
both superposed in one coordinate frame via the same C-terminal
peptide-backbone Kabsch alignment aligned_pair_figure.py uses
(lib/peptide_align.py's `load_and_align_pair`) — but drawn together in one
panel instead of two side-by-side panels, so the "gross change" to the
domain and its peptide interface reads directly as an overlay. Four
independent colors distinguish which structure is which: reference domain /
reference peptide / design domain / design peptide (defaults: cyan / green /
red / yellow).

As in aligned_pair_figure.py, `design` is oriented with the standard
domain/peptide long-axis-vertical convention and `reference` is rigidly
Kabsch-fit so its peptide's C-terminal binding motif (the last
`--align-length` residues) lands on `design`'s equivalent residues — see
that script's docstring for the full rationale. The difference here is
purely in rendering: both objects stay visible together and the camera pans
to fit the whole overlaid pair (never a tighter crop, so nothing overlaid
extends past the canvas edge — same rationale as
motif_superposition_figure.py's "reference or design" pan/zoom target).

An overlay is easy to misread when two cartoons occlude each other, so this
always offers (via optional `--out-reference`/`--out-design`) the same
solo, unobstructed views aligned_pair_figure.py's two-panel figure produces
— every overlay-style figure script should (see lib/pymol_scene.py's
`render_solo`, shared by both scripts).

Usage:
    python3 aligned_overlay_figure.py <reference.pdb> <design.pdb> <output.png> \
        [--reference-chain-domain A] [--reference-chain-peptide B] \
        [--design-chain-domain A] [--design-chain-peptide B] \
        [--align-length 4] \
        [--reference-domain-color cyan] [--reference-peptide-color green] \
        [--design-domain-color red] [--design-peptide-color yellow] \
        [--zoom-buffer 5] \
        [--out-reference reference_alone.png] [--out-design design_alone.png] \
        [--width 1800] [--height 2400] [--dpi 300] [--bg white]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from peptide_align import load_and_align_pair  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save, render_solo  # noqa: E402


def style_overlay(design_obj, design_domain_chain, design_peptide_chain, design_domain_color, design_peptide_color,
                   reference_obj, reference_domain_chain, reference_peptide_chain, reference_domain_color, reference_peptide_color):
    cmd.hide("everything", f"{design_obj} or {reference_obj}")

    cmd.show("cartoon", f"{design_obj} and chain {design_domain_chain}")
    cmd.color(design_domain_color, f"{design_obj} and chain {design_domain_chain}")
    cmd.show("cartoon", f"{design_obj} and chain {design_peptide_chain}")
    cmd.color(design_peptide_color, f"{design_obj} and chain {design_peptide_chain}")

    cmd.show("cartoon", f"{reference_obj} and chain {reference_domain_chain}")
    cmd.color(reference_domain_color, f"{reference_obj} and chain {reference_domain_chain}")
    cmd.show("cartoon", f"{reference_obj} and chain {reference_peptide_chain}")
    cmd.color(reference_peptide_color, f"{reference_obj} and chain {reference_peptide_chain}")


def build_figure(reference_pdb, design_pdb, out_png,
                  reference_chain_domain, reference_chain_peptide,
                  design_chain_domain, design_chain_peptide,
                  align_length,
                  reference_domain_color, reference_peptide_color,
                  design_domain_color, design_peptide_color,
                  zoom_buffer, out_reference_png, out_design_png,
                  width, height, dpi, bg):
    load_and_align_pair(
        reference_pdb, design_pdb, reference_chain_peptide,
        design_chain_domain, design_chain_peptide,
        align_length,
    )

    style_overlay(
        "design", design_chain_domain, design_chain_peptide, design_domain_color, design_peptide_color,
        "reference", reference_chain_domain, reference_chain_peptide, reference_domain_color, reference_peptide_color,
    )
    apply_material_aoshiny()

    # Solo, unobstructed views of each structure -- optional, but always
    # available alongside the overlay (see the module docstring).
    if out_reference_png:
        render_solo("reference", ["design"], out_reference_png, zoom_buffer, width, height, dpi, bg)
    if out_design_png:
        render_solo("design", ["reference"], out_design_png, zoom_buffer, width, height, dpi, bg)

    cmd.bg_color(bg)
    # Pan/zoom to fit both structures together, never a tighter crop, so
    # neither overlaid structure ends up cut off at the canvas edge.
    cmd.zoom("reference or design", buffer=zoom_buffer)

    ray_trace_and_save(out_png, width, height, dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference_pdb", help="Design-template structure (e.g. a prod_in/<base>_in/*.pdb scaffold)")
    parser.add_argument("design_pdb", help="Resulting complex structure")
    parser.add_argument("output_png")
    parser.add_argument("--reference-chain-domain", default="A", help="Domain chain id in reference_pdb (default: A)")
    parser.add_argument("--reference-chain-peptide", default="B", help="Peptide chain id in reference_pdb (default: B)")
    parser.add_argument("--design-chain-domain", default="A", help="Domain chain id in design_pdb (default: A)")
    parser.add_argument("--design-chain-peptide", default="B", help="Peptide chain id in design_pdb (default: B)")
    parser.add_argument("--align-length", type=int, default=4,
                         help="Number of C-terminal peptide residues used for the Kabsch fit that places "
                              "reference_pdb into design_pdb's coordinate frame (default: 4, e.g. a PDZ "
                              "target's conserved 'EPEA' binding motif)")
    parser.add_argument("--reference-domain-color", default="cyan", help="PyMOL color for the reference domain cartoon (default: cyan)")
    parser.add_argument("--reference-peptide-color", default="green", help="PyMOL color for the reference peptide cartoon (default: green)")
    parser.add_argument("--design-domain-color", default="red", help="PyMOL color for the design domain cartoon (default: red)")
    parser.add_argument("--design-peptide-color", default="yellow", help="PyMOL color for the design peptide cartoon (default: yellow)")
    parser.add_argument("--zoom-buffer", type=float, default=5.0, help="Padding (Angstroms) when panning/zooming to fit both structures together (default: 5.0)")
    parser.add_argument("--out-reference", default=None, help="Also write a solo, unobstructed render of just the reference structure to this path (default: not written)")
    parser.add_argument("--out-design", default=None, help="Also write a solo, unobstructed render of just the design structure to this path (default: not written)")
    add_render_flags(parser)
    args = parser.parse_args()

    build_figure(
        args.reference_pdb, args.design_pdb, args.output_png,
        args.reference_chain_domain, args.reference_chain_peptide,
        args.design_chain_domain, args.design_chain_peptide,
        args.align_length,
        args.reference_domain_color, args.reference_peptide_color,
        args.design_domain_color, args.design_peptide_color,
        args.zoom_buffer, args.out_reference, args.out_design,
        args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.output_png}")
    if args.out_reference:
        print(f"Wrote {args.out_reference}")
    if args.out_design:
        print(f"Wrote {args.out_design}")


if __name__ == "__main__":
    main()
