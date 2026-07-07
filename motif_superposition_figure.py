#!/usr/bin/env python3
"""
GENERATE step: render a folded design's protein (cartoon) and its motif hot
spot (spheres), Kabsch-aligned onto a reference active-site structure whose
ligand is drawn as licorice. Illustrates how closely a discontinuous-
scaffold design recapitulates a native active site: the reference here
(see the discontinuous-scaffolds project's `scripts/mcsa_41/*.pdb`) is an
M-CSA-style minimal extract — the handful of fixed motif residues plus the
ligand, not a full chain — so it isn't itself cartoon-renderable. The full,
cartoon-worthy protein is the *design's* folded structure (a Chai-1 output);
its "select_fixed_atoms" motif residues are the ones RFDiffusion3 held
fixed at the reference's coordinates during generation, i.e. the binder's
interface hot spot, shown as spheres once the design is aligned back into
the reference's frame (Chai-1 refolds from sequence, so its raw output
coordinates are in an arbitrary frame of their own).

Reads two things per model:
  - `panel_spec` JSON: {"model_name", "reference_pdb", "design_json",
    "design_cif", ...any extra metadata (motif_rmsd, island_count, ...)}
    (see examples/discontinuous_scaffolds_motif/find_best_fold.py, which
    produces these from a completed discontinuous-scaffolds campaign)
  - the RFD3 structured design-spec JSON named by `design_json`, for the
    model's `contig`, `ligand`, and `select_fixed_atoms` (lib/rfd3_motif_select.py)

An optional `anchor_positions` field in `panel_spec` (a list of chai
sequence positions, e.g. `[36]`) splits the motif hot-spot highlight into two
groups instead of one: motif residues at one of those positions are colored
`--anchor-color` (default burnt orange) instead of `--hotspot-color`, for
highlighting which part of the motif is an anchor residue — one identified
by the discontinuous-scaffolds project's analysis.py as well-predicted
enough to hold fixed in a subsequent redesign, or (since chai sequence
positions are invariant across redesign generations — see
lib/rfd3_motif_select.py's `anchor_chai_positions`) the corresponding
carried-over portion of a later redesign generation's own motif. Omitting
the field (or leaving it empty) reproduces the original single-color
behavior exactly. See examples/discontinuous_scaffolds_anchor_progression/
for a worked two-panel (initial design / final redesign) pipeline built on
this.

Both groups are shown per `--motif-representation`: `spheres` (default,
unchanged from the original behavior) draws each group's named
`select_fixed_atoms` atoms as spheres, colored `--hotspot-color`/
`--anchor-color` respectively; `cartoon` instead just colors each group's
residues' whole cartoon segment, with no added spheres at all — used by
examples/discontinuous_scaffolds_anchor_progression/ for a spheres-free
figure where the anchor highlight is a colored cartoon segment like the
rest of the motif.

The reference structure and the design fold are two independent
coordinate frames — the Kabsch alignment (lib/kabsch.py) that places the
design's motif into the reference's frame is computed from backbone (N,
CA, C, O) atoms of every motif residue, using the contig's chai-sequence-
position mapping to find each residue's equivalent atoms in the design
fold; the actual hot-spot atoms shown (as spheres or a colored cartoon
segment, see above) are whatever atom names `select_fixed_atoms` names for
that residue (often a sidechain subset) or, in `cartoon` mode, the whole
residue.

Usage:
    python3 motif_superposition_figure.py <panel_spec.json> <output.png> \
        [--protein-color yellow] [--hotspot-color purple] [--ligand-color blue] \
        [--motif-representation spheres|cartoon] [--sphere-scale 0.5] \
        [--ligand-representation licorice|surface] \
        [--ligand-transparency 0.5] [--orient-toward ligand|hotspot] \
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
from orient import orient_long_axis_vertical, orient_look_down_on_plane  # noqa: E402
from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save  # noqa: E402


def compute_motif_alignment(protein_atoms, contig_map):
    """Return (R, t) aligning the 'design' object onto the 'reference'
    object, fit from backbone (N, CA, C, O) atoms of every motif residue in
    `protein_atoms` (keys are (chain, resnum) in reference numbering)."""
    ref_pts, design_pts = [], []
    for (chain, resnum) in protein_atoms:
        chai_pos = contig_map.get((chain, resnum))
        if chai_pos is None:
            continue
        for atom_name in ("N", "CA", "C", "O"):
            ref_sel = f"reference and chain {chain} and resi {resnum} and name {atom_name}"
            design_sel = f"design and resi {chai_pos} and name {atom_name}"
            if cmd.count_atoms(ref_sel) != 1 or cmd.count_atoms(design_sel) != 1:
                continue
            ref_pts.append(cmd.get_coords(ref_sel)[0])
            design_pts.append(cmd.get_coords(design_sel)[0])

    if len(ref_pts) < 3:
        sys.exit(f"Only {len(ref_pts)} matched backbone atom pair(s) across the motif residues; need >= 3 for Kabsch alignment")

    return kabsch_fit(np.array(design_pts), np.array(ref_pts))


BACKBONE_ATOMS = ("N", "CA", "C", "O")

# PyMOL has no built-in "burnt orange" (its stock "orange"/"tv_orange" are
# brighter and more saturated); register the standard burnt-orange hex
# (#CC5500) as a named color so --anchor-color burnt_orange (the default)
# works like any built-in color name.
BURNT_ORANGE_RGB = (0.8, 0.333, 0.0)


def build_hotspot_selection(protein_atoms, contig_map, anchor_positions=None, by_residue=False):
    """PyMOL selection(s) (on the 'design' object, post-alignment) for the
    motif residues named in select_fixed_atoms, at their chai sequence
    positions, split into an anchor group (chai sequence position in
    `anchor_positions`) and everything else.

    When `by_residue` is False (the original, sphere-rendered behavior),
    both groups are selected by their named atoms: `select_fixed_atoms`
    names the atoms RFDiffusion3 held fixed during backbone generation, but
    the downstream LigandMPNN sequence design step is free to redesign a
    position's residue identity unless that position was *also* pinned for
    sequence design — so a motif residue named only by its sidechain atoms
    (e.g. an Arg's "NE,CD,CZ,CG,NH1,NH2") can still end up a different
    residue in the folded design. When none of the named atoms survive at a
    residue, fall back to its backbone atoms (always present) so the
    position still reads as a hot spot rather than silently vanishing.

    When `by_residue` is True, both groups are instead selected by *whole
    residue* — the right shape for coloring a cartoon segment rather than
    drawing spheres, since cartoon rendering doesn't care which specific
    atoms are named.

    Returns (anchor_sel, hotspot_sel): motif residues whose chai sequence
    position is in `anchor_positions` go into `anchor_sel`, everything else
    into `hotspot_sel`. With `anchor_positions` empty/None, `anchor_sel` is
    always "none" and every motif residue falls into `hotspot_sel`."""
    anchor_positions = anchor_positions or set()
    anchor_terms, hotspot_terms = [], []
    for (chain, resnum), atom_names in protein_atoms.items():
        chai_pos = contig_map.get((chain, resnum))
        if chai_pos is None:
            continue
        is_anchor = chai_pos in anchor_positions
        if by_residue:
            term = f"(design and resi {chai_pos})"
        else:
            if not atom_names:
                continue
            present = [a for a in atom_names if cmd.count_atoms(f"design and resi {chai_pos} and name {a}") == 1]
            term = f"(design and resi {chai_pos} and name {'+'.join(present or BACKBONE_ATOMS)})"
        (anchor_terms if is_anchor else hotspot_terms).append(term)
    anchor_sel = " or ".join(anchor_terms) if anchor_terms else "none"
    hotspot_sel = " or ".join(hotspot_terms) if hotspot_terms else "none"
    return anchor_sel, hotspot_sel


def show_ligand(ligand_sel, representation, color, transparency):
    """Render the ligand as either bold licorice (the original convention)
    or a transparent surface (for a close-up figure where an opaque or
    stick ligand would hide the motif spheres nested against it)."""
    if representation == "surface":
        # PyMOL silently auto-flags small, disconnected hetero groups (the
        # whole point of the M-CSA-style minimal reference structures this
        # script loads) as "ignore" on load, which zeroes their surface
        # area with no error — clear it or `show surface` renders nothing.
        cmd.flag("ignore", ligand_sel, "clear")
        cmd.set("surface_quality", 2)
        # Without this, back-facing surface/cartoon polygons (unavoidable
        # this close to a motif nested inside a fold) render solid black
        # instead of shaded-through, which reads as a rendering glitch.
        cmd.set("two_sided_lighting", 1)
        cmd.show("surface", ligand_sel)
        cmd.set("transparency", transparency, ligand_sel)
    else:
        cmd.show("sticks", ligand_sel)
        cmd.set("stick_radius", 0.3, ligand_sel)
        cmd.set("stick_quality", 15)
    cmd.color(color, ligand_sel)


def build_figure(panel_spec_path, out_png, protein_color, hotspot_color, anchor_color, ligand_color,
                  motif_representation, sphere_scale, ligand_representation, ligand_transparency,
                  orient_toward, zoom_buffer, trim, trim_pad, width, height, dpi, bg):
    cmd.set_color("burnt_orange", list(BURNT_ORANGE_RGB))

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

    cmd.load(reference_pdb, "reference")
    cmd.remove("reference and hydro")
    cmd.load(design_cif, "design")
    cmd.remove("design and hydro")

    R, t = compute_motif_alignment(protein_atoms, contig_map)
    cmd.transform_selection("design", as_homogeneous(R, t), homogenous=1)

    anchor_sel, hotspot_sel = build_hotspot_selection(
        protein_atoms, contig_map, anchor_positions, by_residue=(motif_representation == "cartoon"),
    )
    if cmd.count_atoms(anchor_sel) == 0 and cmd.count_atoms(hotspot_sel) == 0:
        sys.exit(f"{model_name}: hot-spot selections matched no atoms in the aligned design")

    # `reference` is an M-CSA-style minimal active-site extract (the fixed
    # motif residues + ligand, ~a dozen atoms) rather than a full chain, so
    # it isn't cartoon-renderable; `design` is the actual full-length folded
    # protein, now Kabsch-aligned into the reference's coordinate frame.
    protein_sel = "design and polymer.protein"
    ligand_sel = ligand_resn_selection("reference", ligand_resns)
    if cmd.count_atoms(ligand_sel) == 0:
        sys.exit(f"{model_name}: ligand selection {ligand_sel!r} matched no atoms in {reference_pdb}")

    motif_sel = f"({anchor_sel}) or ({hotspot_sel})"

    cmd.dss("design")
    if orient_toward == "hotspot":
        orient_look_down_on_plane(motif_sel, protein_sel, "reference or design")
    else:
        orient_long_axis_vertical(f"{protein_sel} and name CA and polymer", ligand_sel, "reference or design")

    cmd.hide("everything", "reference or design")

    cmd.show("cartoon", protein_sel)
    cmd.color(protein_color, protein_sel)

    show_ligand(ligand_sel, ligand_representation, ligand_color, ligand_transparency)

    if cmd.count_atoms(hotspot_sel) > 0:
        cmd.color(hotspot_color, hotspot_sel)
        if motif_representation == "spheres":
            cmd.show("spheres", hotspot_sel)
            cmd.set("sphere_scale", sphere_scale, hotspot_sel)

    if cmd.count_atoms(anchor_sel) > 0:
        cmd.color(anchor_color, anchor_sel)
        if motif_representation == "spheres":
            cmd.show("spheres", anchor_sel)
            cmd.set("sphere_scale", sphere_scale, anchor_sel)

    apply_material_aoshiny()

    cmd.bg_color(bg)
    # Pan/zoom to fit everything actually shown (cartoon + spheres + ligand)
    # first, so nothing is cropped out of the 3D render itself — a tight
    # zoom straight to just the hot-spot/ligand region can leave the wider
    # cartoon extending past the canvas edge (seen as the structure being
    # "cut off"). Any tighter, close-up framing is done safely afterward, as
    # a 2D crop of the fully-rendered image (see `trim` below), not by
    # narrowing the 3D camera.
    cmd.zoom("reference or design", buffer=zoom_buffer)

    ray_trace_and_save(out_png, width, height, dpi)

    if trim:
        trimmed = trim_to_content(Image.open(out_png), bg=bg, pad=trim_pad)
        trimmed.save(out_png, dpi=(dpi, dpi))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("panel_spec", help="Path to a panel-spec JSON (model_name, reference_pdb, design_json, design_cif, ...)")
    parser.add_argument("output_png")
    parser.add_argument("--protein-color", default="yellow", help="PyMOL color for the reference protein cartoon (default: yellow)")
    parser.add_argument("--hotspot-color", default="purple", help="PyMOL color for the design's motif hot-spot atom spheres (default: purple)")
    parser.add_argument("--anchor-color", default="burnt_orange",
                         help="PyMOL color for motif atoms named in the panel spec's optional 'anchor_positions' "
                              "field (default: burnt_orange, a custom #CC5500 registered at run time since PyMOL "
                              "has no built-in burnt orange); atoms not in that list keep --hotspot-color")
    parser.add_argument("--ligand-color", default="blue", help="PyMOL color for the reference ligand (default: blue)")
    parser.add_argument("--motif-representation", default="spheres", choices=["spheres", "cartoon"],
                         help="How to highlight every motif residue (both hot-spot and, if the panel spec sets "
                              "'anchor_positions', anchor residues): 'spheres' (default, original behavior) draws "
                              "each group's select_fixed_atoms atoms as spheres, colored --hotspot-color/"
                              "--anchor-color respectively; 'cartoon' instead just colors each group's whole "
                              "cartoon segment, with no spheres at all.")
    parser.add_argument("--sphere-scale", type=float, default=0.5, help="Sphere radius scale for the hot-spot atoms (default: 0.5)")
    parser.add_argument("--ligand-representation", default="licorice", choices=["licorice", "surface"],
                         help="How to render the reference ligand (default: licorice)")
    parser.add_argument("--ligand-transparency", type=float, default=0.0,
                         help="Surface transparency, 0 (opaque) to 1 (invisible); only applies to --ligand-representation surface (default: 0.0)")
    parser.add_argument("--orient-toward", default="ligand", choices=["ligand", "hotspot"],
                         help="How to orient the scene: 'ligand' (default) keeps the protein's long axis vertical "
                              "and twists the ligand to face the camera, for a whole-complex panel; 'hotspot' looks "
                              "straight down onto the best-fit plane of the motif hot-spot atoms, for a close-up "
                              "panel that faces the binding interface itself")
    parser.add_argument("--zoom-buffer", type=float, default=5.0,
                         help="Padding (Angstroms) around the full reference+design complex when panning/zooming "
                              "the camera to fit it — always the whole complex, never just the motif, so nothing "
                              "shown ends up cropped out of the 3D render (default: 5.0)")
    parser.add_argument("--no-trim", action="store_true",
                         help="Skip cropping the rendered PNG down to its content bounding box (by default this "
                              "runs after the pan/zoom above, removing excess background margin for a tighter panel)")
    parser.add_argument("--trim-pad", type=int, default=20, help="Pixels of background left around the content after trimming (default: 20)")
    add_render_flags(parser, default_height=1800)
    args = parser.parse_args()

    build_figure(
        args.panel_spec, args.output_png, args.protein_color, args.hotspot_color, args.anchor_color,
        args.ligand_color, args.motif_representation, args.sphere_scale, args.ligand_representation,
        args.ligand_transparency, args.orient_toward, args.zoom_buffer, not args.no_trim, args.trim_pad,
        args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
