#!/usr/bin/env python3
"""
GENERATE step: render a de novo RFDiffusion3 backbone against the downstream
fold prediction (AlphaFold2/ColabFold) it led to, in production
small-molecule-binder campaigns of the kind under
prod/small_molecule_binding-*/p<N>/<task>_rfd3 -> .../<task>_alphafold.

The fold prediction is typically of the *whole* oligomeric assembly the
backbone seeds, folded as one long single-chain object (N copies of the
protomer, concatenated, continuous residue numbering) rather than the lone
protomer the backbone itself is -- so this uses lib/oligomer_align.py to
find and align onto whichever window of the fold actually corresponds to
the backbone, rather than a plain whole-structure cealign.

Four images are produced:
  - the backbone alone (--alt-representation, default surface)
  - the fold alone (cartoon)
  - both aligned and overlapped in one panel (backbone in --alt-representation,
    fold cartoon)
  - the same overlap with representations/colors swapped (backbone cartoon,
    fold in --alt-representation) -- a second read of the same alignment

`--alt-representation` (default: surface) is whichever non-cartoon PyMOL
representation is being compared against cartoon across a pair of figures;
originally ribbon, now surface, without duplicating the whole script per
representation choice.

Optionally (`--out-fold-alt`), a fifth image renders the fold alone in the
second pass's --alt-representation -- the same one shown alongside the
backbone in the swapped overlap, but by itself.

Two more fold-alone images are available, independent of alt_representation:
  - `--out-fold-ribbon`: the fold alone in ribbon, --fold-color (default red)
  - `--out-fold-ss`: the fold alone in licorice, colored by secondary
    structure (--helix-color/--sheet-color/--loop-color) rather than one
    flat color, so helix/sheet/loop content reads directly off the fold

Usage:
    python3 backbone_fold_overlay_figure.py <backbone.cif> <fold.pdb> \
        <out_backbone.png> <out_fold.png> <out_overlap.png> <out_overlap_swapped.png> \
        [--backbone-chain A] [--fold-chain A] \
        [--backbone-color blue] [--fold-color red] [--swapped-backbone-color green] \
        [--alt-representation surface] \
        [--zoom-buffer 5] [--out-fold-alt out_fold_alt.png] \
        [--out-fold-ribbon out_fold_ribbon.png] [--out-fold-ss out_fold_ss.png] \
        [--helix-color red] [--sheet-color yellow] [--loop-color green] \
        [--width 1800] [--height 2400] [--dpi 300] [--bg white]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from oligomer_align import align_onto_best_repeat  # noqa: E402
from orient import orient_long_axis_vertical  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save, render_solo  # noqa: E402


def load_and_align(backbone_path, fold_path, backbone_chain, fold_chain,
                    backbone_obj="design_backbone", fold_obj="predicted_fold"):
    cmd.load(backbone_path, backbone_obj)
    cmd.remove(f"{backbone_obj} and hydro")
    cmd.remove(f"{backbone_obj} and not chain {backbone_chain}")
    cmd.dss(backbone_obj)

    cmd.load(fold_path, fold_obj)
    cmd.remove(f"{fold_obj} and hydro")
    cmd.dss(fold_obj)

    backbone_ca_sel = f"{backbone_obj} and polymer.protein and name CA"
    window_size = cmd.count_atoms(backbone_ca_sel)
    lo, hi, result = align_onto_best_repeat(backbone_ca_sel, fold_obj, fold_chain, window_size)
    print(
        f"Aligned {backbone_obj} onto {fold_obj}/chain {fold_chain} resi {lo}-{hi}: "
        f"RMSD={result['RMSD']:.3f} over {result['alignment_length']} residues",
        file=sys.stderr,
    )

    # Orient the whole aligned pair together so the fold assembly's own long
    # axis runs vertical -- it's the larger, more complete structure of the
    # two, so it makes the more stable reference for framing both renders.
    orient_long_axis_vertical(
        f"{fold_obj} and chain {fold_chain} and name CA and polymer",
        None,
        f"{backbone_obj} or {fold_obj}",
    )
    return backbone_obj, fold_obj


def style(backbone_obj, backbone_chain, backbone_repr, backbone_color,
          fold_obj, fold_chain, fold_repr, fold_color):
    cmd.hide("everything", f"{backbone_obj} or {fold_obj}")
    cmd.show(backbone_repr, f"{backbone_obj} and chain {backbone_chain}")
    cmd.color(backbone_color, f"{backbone_obj} and chain {backbone_chain}")
    cmd.show(fold_repr, f"{fold_obj} and chain {fold_chain}")
    cmd.color(fold_color, f"{fold_obj} and chain {fold_chain}")


def style_by_secondary_structure(obj, chain, repr_, helix_color, sheet_color, loop_color):
    """Hide/show+color `obj`'s `chain` in `repr_`, split by cmd.dss's own
    per-residue ss assignment (already run in load_and_align) rather than
    one flat color -- helix/sheet/loop content reads directly off the
    figure instead of needing a separate cartoon panel to show it.

    cmd.dss only tags each residue's CA (not its side-chain atoms) with the
    ss code, so a plain atom-level "and ss H" selection would leave every
    side chain uncolored (falling through to loop_color) in an all-atom
    representation like licorice -- `byres` spreads each residue's CA-level
    assignment across all of that residue's atoms first."""
    sel = f"{obj} and chain {chain}"
    cmd.hide("everything", obj)
    cmd.show(repr_, sel)
    cmd.color(helix_color, f"byres ({sel} and ss H)")
    cmd.color(sheet_color, f"byres ({sel} and ss S)")
    cmd.color(loop_color, f"({sel}) and not (byres ({sel} and (ss H or ss S)))")


def build_figure(backbone_path, fold_path, out_backbone_png, out_fold_png,
                  out_overlap_png, out_overlap_swapped_png,
                  backbone_chain, fold_chain,
                  backbone_color, fold_color, swapped_backbone_color,
                  zoom_buffer, width, height, dpi, bg,
                  alt_representation="surface", out_fold_alt_png=None,
                  out_fold_ribbon_png=None, out_fold_ss_png=None,
                  helix_color="red", sheet_color="yellow", loop_color="green"):
    backbone_obj, fold_obj = load_and_align(backbone_path, fold_path, backbone_chain, fold_chain)

    # Pass 1: backbone in alt_representation, fold as cartoon -- used for
    # the two solo views and the first overlap.
    style(backbone_obj, backbone_chain, alt_representation, backbone_color,
          fold_obj, fold_chain, "cartoon", fold_color)
    apply_material_aoshiny()

    render_solo(backbone_obj, [fold_obj], out_backbone_png, zoom_buffer, width, height, dpi, bg)
    render_solo(fold_obj, [backbone_obj], out_fold_png, zoom_buffer, width, height, dpi, bg)

    cmd.bg_color(bg)
    cmd.viewport(width, height)
    cmd.zoom(f"{backbone_obj} or {fold_obj}", buffer=zoom_buffer)
    ray_trace_and_save(out_overlap_png, width, height, dpi)

    # Pass 2: same overlap, backbone as cartoon / fold in alt_representation,
    # backbone recolored -- same camera (zoom target is the same pair of
    # objects).
    style(backbone_obj, backbone_chain, "cartoon", swapped_backbone_color,
          fold_obj, fold_chain, alt_representation, fold_color)
    cmd.viewport(width, height)
    cmd.zoom(f"{backbone_obj} or {fold_obj}", buffer=zoom_buffer)
    ray_trace_and_save(out_overlap_swapped_png, width, height, dpi)

    if out_fold_alt_png:
        render_solo(fold_obj, [backbone_obj], out_fold_alt_png, zoom_buffer, width, height, dpi, bg)

    if out_fold_ribbon_png:
        cmd.hide("everything", f"{backbone_obj} or {fold_obj}")
        cmd.show("ribbon", f"{fold_obj} and chain {fold_chain}")
        cmd.color(fold_color, f"{fold_obj} and chain {fold_chain}")
        render_solo(fold_obj, [backbone_obj], out_fold_ribbon_png, zoom_buffer, width, height, dpi, bg)

    if out_fold_ss_png:
        style_by_secondary_structure(fold_obj, fold_chain, "licorice", helix_color, sheet_color, loop_color)
        render_solo(fold_obj, [backbone_obj], out_fold_ss_png, zoom_buffer, width, height, dpi, bg)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("backbone_path", help="RFDiffusion3 de novo backbone output (e.g. <N>_rfd3/out/*_model_0.cif.gz)")
    parser.add_argument("fold_path", help="Downstream AlphaFold2/ColabFold prediction (e.g. <M>_alphafold/out/*rank_001*.pdb)")
    parser.add_argument("out_backbone_png")
    parser.add_argument("out_fold_png")
    parser.add_argument("out_overlap_png")
    parser.add_argument("out_overlap_swapped_png")
    parser.add_argument("--backbone-chain", default="A", help="Protein chain id in backbone_path (default: A)")
    parser.add_argument("--fold-chain", default="A", help="Protein chain id in fold_path (default: A)")
    parser.add_argument("--backbone-color", default="blue", help="PyMOL color for the backbone in the solo view and first overlap (default: blue)")
    parser.add_argument("--fold-color", default="red", help="PyMOL color for the fold, both overlaps (default: red)")
    parser.add_argument("--swapped-backbone-color", default="green", help="PyMOL color for the backbone in the second (swapped-representation) overlap (default: green)")
    parser.add_argument("--zoom-buffer", type=float, default=5.0, help="Padding (Angstroms) when panning/zooming (default: 5.0)")
    parser.add_argument("--alt-representation", default="surface", help="Representation shown against cartoon in the solo backbone view, first overlap (on the backbone), and swapped overlap (on the fold) (default: surface)")
    parser.add_argument("--out-fold-alt", default=None, help="Also write a solo render of just the fold, in --alt-representation (default: not written)")
    parser.add_argument("--out-fold-ribbon", default=None, help="Also write a solo render of just the fold, in ribbon/--fold-color (default: not written)")
    parser.add_argument("--out-fold-ss", default=None, help="Also write a solo render of just the fold, in licorice colored by secondary structure (default: not written)")
    parser.add_argument("--helix-color", default="red", help="Color for helix residues in --out-fold-ss (default: red)")
    parser.add_argument("--sheet-color", default="yellow", help="Color for sheet residues in --out-fold-ss (default: yellow)")
    parser.add_argument("--loop-color", default="green", help="Color for loop/coil residues in --out-fold-ss (default: green)")
    add_render_flags(parser)
    args = parser.parse_args()

    build_figure(
        args.backbone_path, args.fold_path,
        args.out_backbone_png, args.out_fold_png, args.out_overlap_png, args.out_overlap_swapped_png,
        args.backbone_chain, args.fold_chain,
        args.backbone_color, args.fold_color, args.swapped_backbone_color,
        args.zoom_buffer, args.width, args.height, args.dpi, args.bg,
        alt_representation=args.alt_representation, out_fold_alt_png=args.out_fold_alt,
        out_fold_ribbon_png=args.out_fold_ribbon, out_fold_ss_png=args.out_fold_ss,
        helix_color=args.helix_color, sheet_color=args.sheet_color, loop_color=args.loop_color,
    )
    for path in (args.out_backbone_png, args.out_fold_png, args.out_overlap_png, args.out_overlap_swapped_png):
        print(f"Wrote {path}")
    for path in (args.out_fold_alt, args.out_fold_ribbon, args.out_fold_ss):
        if path:
            print(f"Wrote {path}")


if __name__ == "__main__":
    main()
