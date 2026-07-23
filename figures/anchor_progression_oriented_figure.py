#!/usr/bin/env python3
"""
GENERATE step: render a two-panel "root design vs. redesign" anchor-
progression figure sharing one fixed camera orientation (the mechanism from
anchor_progression_zoned_figure.py), but using motif_superposition_figure.py's
original representation - a single flat-colored protein cartoon plus colored
cartoon segments for hot-spot/anchor motif residues - rather than
anchor_progression_zoned_figure.py's distance-from-motif surface/cartoon
zoning.

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

Coloring: the whole protein cartoon is `--protein-color` (default yellow);
motif hot-spot residues are `--hotspot-color` (default pink); anchor
residues (the panel spec's optional `anchor_positions` field) are
`--anchor-color` (default green).

Reads two panel-spec JSONs (see motif_superposition_figure.py's docstring
for the schema; both produced by find_anchor_progression.py): `initial_spec`
(the root design - also the source of the shared camera orientation) and
`final_spec` (the passing redesign).

Usage:
    python3 anchor_progression_oriented_figure.py <initial_spec.json> <final_spec.json> \
        <out_initial.png> <out_final.png> \
        [--protein-color yellow] [--hotspot-color pink] [--anchor-color green] \
        [--ligand-color blue] [--ligand-representation licorice|surface] [--ligand-transparency 0.0] \
        [--zoom-buffer 5.0] [--no-trim] [--trim-pad 20] \
        [--width 1800] [--height 1800] [--dpi 300] [--bg white]
"""
import argparse
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
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
    motif is always shown as a cartoon segment here (never spheres), so -
    unlike motif_superposition_figure.py's build_hotspot_selection - there's
    no named-atom/spheres fallback to handle."""
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


def render_panel(panel, anchor_sel, hotspot_sel, out_png,
                  protein_color, hotspot_color, anchor_color, ligand_color,
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

    cmd.show("cartoon", protein_sel)
    cmd.color(protein_color, protein_sel)

    show_ligand(ligand_sel, ligand_representation, ligand_color, ligand_transparency)

    if cmd.count_atoms(hotspot_sel) > 0:
        cmd.color(hotspot_color, hotspot_sel)

    if cmd.count_atoms(anchor_sel) > 0:
        cmd.color(anchor_color, anchor_sel)

    apply_material_aoshiny()
    cmd.bg_color(bg)
    cmd.zoom(f"{reference_obj} or {design_obj}", buffer=zoom_buffer)
    ray_trace_and_save(out_png, width, height, dpi)

    if trim:
        trimmed = trim_to_content(Image.open(out_png), bg=bg, pad=trim_pad)
        trimmed.save(out_png, dpi=(dpi, dpi))


def build_figure(initial_spec_path, final_spec_path, out_initial_png, out_final_png,
                  protein_color, hotspot_color, anchor_color, ligand_color,
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

    render_panel(init, init_anchor_sel, init_hotspot_sel, out_initial_png,
                 protein_color, hotspot_color, anchor_color, ligand_color,
                 ligand_representation, ligand_transparency,
                 zoom_buffer, trim, trim_pad, width, height, dpi, bg)

    render_panel(final, final_anchor_sel, final_hotspot_sel, out_final_png,
                 protein_color, hotspot_color, anchor_color, ligand_color,
                 ligand_representation, ligand_transparency,
                 zoom_buffer, trim, trim_pad, width, height, dpi, bg)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("initial_spec", help="Path to the initial (root) panel-spec JSON - also the source of the shared camera orientation")
    parser.add_argument("final_spec", help="Path to the final (passing redesign) panel-spec JSON")
    parser.add_argument("out_initial_png")
    parser.add_argument("out_final_png")
    parser.add_argument("--protein-color", default="yellow", help="PyMOL color for the whole protein cartoon (default: yellow)")
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
        args.protein_color, args.hotspot_color, args.anchor_color, args.ligand_color,
        args.ligand_representation, args.ligand_transparency,
        args.zoom_buffer, not args.no_trim, args.trim_pad,
        args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.out_initial_png}")
    print(f"Wrote {args.out_final_png}")


if __name__ == "__main__":
    main()
