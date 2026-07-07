#!/usr/bin/env python3
"""
GENERATE step: render a two-panel "root design vs. redesign" anchor-
progression figure sharing one fixed camera orientation and a
distance-from-motif coloring scheme, instead of
motif_superposition_figure.py's independently-oriented, flat-protein-color
panels.

Both panel specs are loaded into one PyMOL session so they can share a
*single* rigid camera orientation: the best-fit plane through the *initial*
(root) design's own motif residues is computed once (see
lib/orient.py's compute_plane_orientation, the "look down on the motif
plane" convention), and the exact same (rotation, translation) transform is
then reapplied verbatim to the final panel's reference/design objects - not
independently re-fit from the final design's own (differently-folded) motif
atoms - so a reader's eye stays anchored to one consistent camera framing
across both panels. This works because orientation here is *baked into atom
coordinates* via cmd.transform_selection rather than being a camera move
(same convention aligned_pair_figure.py documents), and because both panel
specs share the same reference_pdb file (see find_anchor_progression.py):
reapplying the identical transform to the final panel's reference object
lands it exactly where the initial panel's reference ended up, since both
start from identical un-transformed coordinates.

Coloring is by distance from the motif rather than a single flat protein
color: protein residues within `--near-cutoff` (default 3 Angstroms) of any
motif residue are shown as opaque cartoon (`--near-color`, default yellow);
everything else is shown as a translucent surface (`--far-color`, default
white, at `--far-alpha` transparency). Motif residues themselves are always
shown as cartoon, colored `--hotspot-color` (default pink), or
`--anchor-color` (default green) for anchor residues - the panel spec's
optional `anchor_positions` field, same convention as
motif_superposition_figure.py.

Reads two panel-spec JSONs (see motif_superposition_figure.py's docstring
for the schema; both produced by find_anchor_progression.py): `initial_spec`
(the root design - also the source of the shared camera orientation) and
`final_spec` (the passing redesign).

Usage:
    python3 anchor_progression_zoned_figure.py <initial_spec.json> <final_spec.json> \
        <out_initial.png> <out_final.png> \
        [--near-color yellow] [--far-color white] [--far-alpha 0.3] [--near-cutoff 10.0] \
        [--hotspot-color pink] [--anchor-color green] [--ligand-color blue] \
        [--ligand-representation licorice|surface] [--ligand-transparency 0.0] \
        [--zoom-buffer 5.0] [--no-trim] [--trim-pad 20] \
        [--width 1800] [--height 1800] [--dpi 300] [--bg white]
"""
import argparse
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from kabsch import kabsch_fit, as_homogeneous  # noqa: E402
from rfd3_motif_select import (  # noqa: E402
    load_design_spec, parse_atom_names, split_fixed_atoms,
    contig_to_chai_positions, ligand_resn_selection,
)
from imgtrim import trim_to_content  # noqa: E402
from orient import compute_plane_orientation  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save  # noqa: E402

BACKBONE_ATOMS = ("N", "CA", "C", "O")


def load_panel(panel_spec_path, suffix):
    """Load a panel spec's reference+design structures as `reference_<suffix>`/
    `design_<suffix>` PyMOL objects (distinct per-panel names, since both
    panels' objects live in the same session at once - see module
    docstring), returning everything else needed to align/select/render it."""
    with open(panel_spec_path) as f:
        panel = json.load(f)

    model_name = panel["model_name"]
    entry = load_design_spec(panel["design_json"], model_name)
    reference_pdb = panel.get("reference_pdb") or entry["input"]
    design_cif = panel["design_cif"]
    anchor_positions = set(panel.get("anchor_positions") or [])

    for path in (reference_pdb, design_cif):
        if not os.path.exists(path):
            sys.exit(f"{model_name}: structure file not found: {path}")

    ligand_resns = parse_atom_names(entry["ligand"])
    contig_map = contig_to_chai_positions(entry["contig"])
    protein_atoms, _ligand_atoms = split_fixed_atoms(entry["select_fixed_atoms"], ligand_resns)

    reference_obj, design_obj = f"reference_{suffix}", f"design_{suffix}"
    cmd.load(reference_pdb, reference_obj)
    cmd.remove(f"{reference_obj} and hydro")
    cmd.load(design_cif, design_obj)
    cmd.remove(f"{design_obj} and hydro")

    return {
        "model_name": model_name, "reference_obj": reference_obj, "design_obj": design_obj,
        "ligand_resns": ligand_resns, "contig_map": contig_map, "protein_atoms": protein_atoms,
        "anchor_positions": anchor_positions,
    }


def compute_motif_alignment(reference_obj, design_obj, protein_atoms, contig_map):
    """Return (R, t) aligning `design_obj` onto `reference_obj`, fit from
    backbone (N, CA, C, O) atoms of every motif residue (see
    motif_superposition_figure.py's function of the same name)."""
    ref_pts, design_pts = [], []
    for (chain, resnum) in protein_atoms:
        chai_pos = contig_map.get((chain, resnum))
        if chai_pos is None:
            continue
        for atom_name in BACKBONE_ATOMS:
            ref_sel = f"{reference_obj} and chain {chain} and resi {resnum} and name {atom_name}"
            design_sel = f"{design_obj} and resi {chai_pos} and name {atom_name}"
            if cmd.count_atoms(ref_sel) != 1 or cmd.count_atoms(design_sel) != 1:
                continue
            ref_pts.append(cmd.get_coords(ref_sel)[0])
            design_pts.append(cmd.get_coords(design_sel)[0])

    if len(ref_pts) < 3:
        sys.exit(f"{design_obj}: only {len(ref_pts)} matched backbone atom pair(s) across the motif residues; need >= 3 for Kabsch alignment")

    return kabsch_fit(np.array(design_pts), np.array(ref_pts))


def build_motif_selections(design_obj, protein_atoms, contig_map, anchor_positions):
    """Whole-residue (design_obj) selections (anchor_sel, hotspot_sel): the
    motif is always shown as cartoon here (never spheres), so - unlike
    motif_superposition_figure.py's build_hotspot_selection - there's no
    named-atom/spheres fallback to handle."""
    anchor_terms, hotspot_terms = [], []
    for (chain, resnum) in protein_atoms:
        chai_pos = contig_map.get((chain, resnum))
        if chai_pos is None:
            continue
        term = f"({design_obj} and resi {chai_pos})"
        (anchor_terms if chai_pos in anchor_positions else hotspot_terms).append(term)
    anchor_sel = " or ".join(anchor_terms) if anchor_terms else "none"
    hotspot_sel = " or ".join(hotspot_terms) if hotspot_terms else "none"
    return anchor_sel, hotspot_sel


def show_ligand(ligand_sel, representation, color, transparency):
    """Render the ligand as either bold licorice or a transparent surface
    (see motif_superposition_figure.py's function of the same name)."""
    if representation == "surface":
        cmd.flag("ignore", ligand_sel, "clear")
        cmd.set("surface_quality", 2)
        cmd.set("two_sided_lighting", 1)
        cmd.show("surface", ligand_sel)
        cmd.set("transparency", transparency, ligand_sel)
    else:
        cmd.show("sticks", ligand_sel)
        cmd.set("stick_radius", 0.3, ligand_sel)
        cmd.set("stick_quality", 15)
    cmd.color(color, ligand_sel)


def render_panel(panel, motif_sel, anchor_sel, hotspot_sel, out_png,
                  near_color, far_color, far_alpha, near_cutoff,
                  hotspot_color, anchor_color, ligand_color,
                  ligand_representation, ligand_transparency,
                  zoom_buffer, trim, trim_pad, width, height, dpi, bg):
    design_obj, reference_obj = panel["design_obj"], panel["reference_obj"]
    protein_sel = f"{design_obj} and polymer.protein"
    ligand_sel = ligand_resn_selection(reference_obj, panel["ligand_resns"])
    if cmd.count_atoms(ligand_sel) == 0:
        sys.exit(f"{panel['model_name']}: ligand selection {ligand_sel!r} matched no atoms")

    # Every other object in the shared session must stay hidden while this
    # panel renders (see module docstring: both panels' objects live in one
    # PyMOL session at once, for shared-orientation purposes).
    cmd.hide("everything", "all")

    near_sel = f"byres (({protein_sel}) within {near_cutoff} of ({motif_sel}))"
    # near_sel is already guaranteed to cover every motif/anchor residue (each
    # is trivially within 0 of itself), but excluding motif_sel here too,
    # explicitly, means far_sel can never include a motif/anchor residue even
    # if near_sel's definition changes later - the surface is meant only as
    # a reference to the protein's background global structure, never the
    # motif itself.
    far_sel = f"({protein_sel}) and not ({near_sel}) and not ({motif_sel})"

    cmd.show("cartoon", near_sel)
    cmd.color(near_color, near_sel)

    if cmd.count_atoms(far_sel) > 0:
        # surface_quality/two_sided_lighting are object-state-level settings
        # (PyMOL warns and silently ignores the selection scope if given a
        # third argument here - see show_ligand's identical surface-mode
        # settings, which are also set globally, not per-selection); without
        # two_sided_lighting in particular, the torn edge where this partial
        # surface meets the unsurfaced near_sel cartoon renders as solid
        # black instead of shaded-through.
        cmd.set("surface_quality", 2)
        cmd.set("two_sided_lighting", 1)
        cmd.show("surface", far_sel)
        cmd.set("transparency", far_alpha, far_sel)
        cmd.color(far_color, far_sel)

    show_ligand(ligand_sel, ligand_representation, ligand_color, ligand_transparency)

    if cmd.count_atoms(hotspot_sel) > 0:
        cmd.show("cartoon", hotspot_sel)
        cmd.color(hotspot_color, hotspot_sel)

    if cmd.count_atoms(anchor_sel) > 0:
        cmd.show("cartoon", anchor_sel)
        cmd.color(anchor_color, anchor_sel)

    apply_material_aoshiny()
    cmd.bg_color(bg)
    cmd.zoom(f"{reference_obj} or {design_obj}", buffer=zoom_buffer)
    ray_trace_and_save(out_png, width, height, dpi)

    if trim:
        trimmed = trim_to_content(Image.open(out_png), bg=bg, pad=trim_pad)
        trimmed.save(out_png, dpi=(dpi, dpi))


def build_figure(initial_spec_path, final_spec_path, out_initial_png, out_final_png,
                  near_color, far_color, far_alpha, near_cutoff,
                  hotspot_color, anchor_color, ligand_color,
                  ligand_representation, ligand_transparency,
                  zoom_buffer, trim, trim_pad, width, height, dpi, bg):
    init = load_panel(initial_spec_path, "init")
    final = load_panel(final_spec_path, "final")

    R, t = compute_motif_alignment(init["reference_obj"], init["design_obj"], init["protein_atoms"], init["contig_map"])
    cmd.transform_selection(init["design_obj"], as_homogeneous(R, t), homogenous=1)

    R, t = compute_motif_alignment(final["reference_obj"], final["design_obj"], final["protein_atoms"], final["contig_map"])
    cmd.transform_selection(final["design_obj"], as_homogeneous(R, t), homogenous=1)

    init_anchor_sel, init_hotspot_sel = build_motif_selections(
        init["design_obj"], init["protein_atoms"], init["contig_map"], init["anchor_positions"])
    final_anchor_sel, final_hotspot_sel = build_motif_selections(
        final["design_obj"], final["protein_atoms"], final["contig_map"], final["anchor_positions"])

    for label, anchor_sel, hotspot_sel in (
        ("initial", init_anchor_sel, init_hotspot_sel), ("final", final_anchor_sel, final_hotspot_sel),
    ):
        if cmd.count_atoms(anchor_sel) == 0 and cmd.count_atoms(hotspot_sel) == 0:
            sys.exit(f"{label}: motif selections matched no atoms in the aligned design")

    init_motif_sel = f"({init_anchor_sel}) or ({init_hotspot_sel})"
    final_motif_sel = f"({final_anchor_sel}) or ({final_hotspot_sel})"

    cmd.dss(f"{init['design_obj']} or {final['design_obj']}")

    # Shared camera: fit the best-fit plane through the *initial* design's
    # own motif only, then reapply the exact same rigid transform to the
    # final panel's objects too - not a fresh fit from its own motif - so
    # both panels share one consistent orientation (see module docstring).
    protein_sel_init = f"{init['design_obj']} and polymer.protein"
    rot, t = compute_plane_orientation(init_motif_sel, protein_sel_init)
    orient_transform = as_homogeneous(rot, t)
    cmd.transform_selection(f"{init['reference_obj']} or {init['design_obj']}", orient_transform, homogenous=1)
    cmd.transform_selection(f"{final['reference_obj']} or {final['design_obj']}", orient_transform, homogenous=1)

    render_panel(init, init_motif_sel, init_anchor_sel, init_hotspot_sel, out_initial_png,
                 near_color, far_color, far_alpha, near_cutoff,
                 hotspot_color, anchor_color, ligand_color,
                 ligand_representation, ligand_transparency,
                 zoom_buffer, trim, trim_pad, width, height, dpi, bg)

    render_panel(final, final_motif_sel, final_anchor_sel, final_hotspot_sel, out_final_png,
                 near_color, far_color, far_alpha, near_cutoff,
                 hotspot_color, anchor_color, ligand_color,
                 ligand_representation, ligand_transparency,
                 zoom_buffer, trim, trim_pad, width, height, dpi, bg)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("initial_spec", help="Path to the initial (root) panel-spec JSON - also the source of the shared camera orientation")
    parser.add_argument("final_spec", help="Path to the final (passing redesign) panel-spec JSON")
    parser.add_argument("out_initial_png")
    parser.add_argument("out_final_png")
    parser.add_argument("--near-color", default="yellow",
                         help="Cartoon color for protein residues within --near-cutoff of the motif (default: yellow)")
    parser.add_argument("--far-color", default="white",
                         help="Surface color for protein residues farther than --near-cutoff from the motif (default: white)")
    parser.add_argument("--far-alpha", type=float, default=0.3,
                         help="Surface transparency for far residues, 0 (opaque) to 1 (invisible) (default: 0.3)")
    parser.add_argument("--near-cutoff", type=float, default=10.0,
                         help="Distance (Angstroms) from the motif within which protein residues are shown near/opaque cartoon rather than far/translucent surface (default: 10.0)")
    parser.add_argument("--hotspot-color", default="pink", help="Cartoon color for motif hot-spot residues (default: pink)")
    parser.add_argument("--anchor-color", default="green",
                         help="Cartoon color for motif residues named in the panel spec's optional 'anchor_positions' field (default: green)")
    parser.add_argument("--ligand-color", default="blue", help="PyMOL color for the reference ligand (default: blue)")
    parser.add_argument("--ligand-representation", default="licorice", choices=["licorice", "surface"],
                         help="How to render the reference ligand (default: licorice)")
    parser.add_argument("--ligand-transparency", type=float, default=0.0,
                         help="Surface transparency, 0 (opaque) to 1 (invisible); only applies to --ligand-representation surface (default: 0.0)")
    parser.add_argument("--zoom-buffer", type=float, default=5.0,
                         help="Padding (Angstroms) when panning/zooming each panel's camera to fit its own reference+design complex (default: 5.0)")
    parser.add_argument("--no-trim", action="store_true",
                         help="Skip cropping each rendered PNG down to its content bounding box")
    parser.add_argument("--trim-pad", type=int, default=20, help="Pixels of background left around the content after trimming (default: 20)")
    add_render_flags(parser, default_height=1800)
    args = parser.parse_args()

    build_figure(
        args.initial_spec, args.final_spec, args.out_initial_png, args.out_final_png,
        args.near_color, args.far_color, args.far_alpha, args.near_cutoff,
        args.hotspot_color, args.anchor_color, args.ligand_color,
        args.ligand_representation, args.ligand_transparency,
        args.zoom_buffer, not args.no_trim, args.trim_pad,
        args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.out_initial_png}")
    print(f"Wrote {args.out_final_png}")


if __name__ == "__main__":
    main()
