#!/usr/bin/env python3
"""
GENERATE step: render a two-panel "template vs. design" comparison row for a
domain/peptide complex — the structure used to seed a design (`reference`,
e.g. a production campaign's prod_in/<model>_in/*.pdb scaffold) and the
resulting complex (`design`, e.g. that campaign's best-scoring predicted
model) — each rendered as its own standalone panel, unlike
motif_superposition_figure.py's single overlaid panel, but sharing one
coordinate frame so a downstream montage of the two panels reads as "the
same peptide, same camera, two different structures around it."

`design` is oriented with the standard domain/peptide convention
(orient_long_axis_vertical: domain long axis vertical, peptide twisted
toward the camera — same as pdz_figure.py/generate_figure.py). `reference`
is *not* independently oriented; instead its coordinates are rigidly
Kabsch-fit (lib/peptide_align.py, this repo's fourth selection/alignment
convention) so its peptide's C-terminal binding motif (the last
`--align-length` residues — e.g. a PDZ target's conserved "EPEA" motif)
lands on `design`'s already-oriented peptide's equivalent C-terminal
residues, matched by position (N -> C order) rather than residue identity —
`reference`'s peptide can be a different total length than `design`'s (a
minimal 4-residue motif vs. a longer full peptide), as long as both share
the same C-terminal motif. Because orientation is baked into atom
coordinates rather than being a camera move (see lib/orient.py), and both
panels are ray-traced from PyMOL's untouched default view, this makes the
shared peptide motif land in the same place on both rendered panels even
though the rest of `reference` (domain sequence/length, extra N-terminal
peptide residues) differs freely from `design`.

Usage:
    python3 aligned_pair_figure.py <reference.pdb> <design.pdb> \
        <out_reference.png> <out_design.png> \
        [--reference-chain-domain A] [--reference-chain-peptide B] \
        [--design-chain-domain A] [--design-chain-peptide B] \
        [--align-length 4] \
        [--domain-color cyan] [--peptide-color green] \
        [--zoom-buffer 5] \
        [--width 1800] [--height 2400] [--dpi 300] [--bg white]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from peptide_align import load_and_align_pair  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, render_solo  # noqa: E402


def style_complex(obj, domain_chain, peptide_chain, domain_color, peptide_color):
    cmd.hide("everything", obj)
    cmd.show("cartoon", f"{obj} and chain {domain_chain}")
    cmd.color(domain_color, f"{obj} and chain {domain_chain}")
    cmd.show("cartoon", f"{obj} and chain {peptide_chain}")
    cmd.color(peptide_color, f"{obj} and chain {peptide_chain}")


def build_figure(reference_pdb, design_pdb, out_reference_png, out_design_png,
                  reference_chain_domain, reference_chain_peptide,
                  design_chain_domain, design_chain_peptide,
                  align_length, domain_color, peptide_color, zoom_buffer,
                  width, height, dpi, bg):
    load_and_align_pair(
        reference_pdb, design_pdb, reference_chain_peptide,
        design_chain_domain, design_chain_peptide,
        align_length,
    )

    style_complex("design", design_chain_domain, design_chain_peptide, domain_color, peptide_color)
    style_complex("reference", reference_chain_domain, reference_chain_peptide, domain_color, peptide_color)
    apply_material_aoshiny()

    render_solo("design", ["reference"], out_design_png, zoom_buffer, width, height, dpi, bg)
    render_solo("reference", ["design"], out_reference_png, zoom_buffer, width, height, dpi, bg)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference_pdb", help="Design-template structure (e.g. a prod_in/<model>_in/*.pdb scaffold), rendered as the reference panel")
    parser.add_argument("design_pdb", help="Resulting complex structure, rendered as the design panel")
    parser.add_argument("out_reference_png")
    parser.add_argument("out_design_png")
    parser.add_argument("--reference-chain-domain", default="A", help="Domain chain id in reference_pdb (default: A)")
    parser.add_argument("--reference-chain-peptide", default="B", help="Peptide chain id in reference_pdb (default: B)")
    parser.add_argument("--design-chain-domain", default="A", help="Domain chain id in design_pdb (default: A)")
    parser.add_argument("--design-chain-peptide", default="B", help="Peptide chain id in design_pdb (default: B)")
    parser.add_argument("--align-length", type=int, default=4,
                         help="Number of C-terminal peptide residues used for the Kabsch fit that places "
                              "reference_pdb into design_pdb's coordinate frame (default: 4, e.g. a PDZ "
                              "target's conserved 'EPEA' binding motif)")
    parser.add_argument("--domain-color", default="cyan", help="PyMOL color for the domain cartoon in both panels (default: cyan)")
    parser.add_argument("--peptide-color", default="green", help="PyMOL color for the peptide cartoon in both panels (default: green)")
    parser.add_argument("--zoom-buffer", type=float, default=5.0, help="Padding (Angstroms) when panning/zooming to fit each panel's own complex (default: 5.0)")
    add_render_flags(parser)
    args = parser.parse_args()

    build_figure(
        args.reference_pdb, args.design_pdb, args.out_reference_png, args.out_design_png,
        args.reference_chain_domain, args.reference_chain_peptide,
        args.design_chain_domain, args.design_chain_peptide,
        args.align_length, args.domain_color, args.peptide_color, args.zoom_buffer,
        args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.out_reference_png}")
    print(f"Wrote {args.out_design_png}")


if __name__ == "__main__":
    main()
