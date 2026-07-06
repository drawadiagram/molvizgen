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

The reference structure and the design fold are two independent
coordinate frames — the Kabsch alignment (lib/kabsch.py) that places the
design's motif into the reference's frame is computed from backbone (N,
CA, C, O) atoms of every motif residue, using the contig's chai-sequence-
position mapping to find each residue's equivalent atoms in the design
fold; the actual "hot spot" atoms shown as spheres are whatever atom names
`select_fixed_atoms` names for that residue (often a sidechain subset).

Usage:
    python3 motif_superposition_figure.py <panel_spec.json> <output.png> \
        [--protein-color yellow] [--hotspot-color purple] [--ligand-color blue] \
        [--sphere-scale 0.5] [--ligand-representation licorice|surface] \
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

import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd  # noqa: E402


def rotation_to_align(a, b):
    """Rotation matrix (3x3) that rotates unit vector a onto unit vector b,
    via Rodrigues' rotation formula. (Shared convention with pdz_figure.py /
    generate_figure.py's orientation step.)"""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        ortho = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, ortho)
        axis /= np.linalg.norm(axis)
        return 2 * np.outer(axis, axis) - np.eye(3)
    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0],
    ])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s ** 2))


def orient_scene_vertically(protein_sel, face_sel, whole_sel):
    """Rotate `whole_sel` rigidly so `protein_sel`'s long axis (largest PCA
    component of its CA atoms) lies vertical (Y), twisted about that axis so
    `face_sel`'s centroid faces the camera (+Z) — the pdz_figure.py /
    generate_figure.py domain/peptide convention, generalized to whatever
    should face the camera: the ligand (the original convention, used when
    the panel is meant to read as a whole complex) or the motif hot-spot
    selection (to bring the binding motif itself into view, e.g. for a
    close-up single-panel figure)."""
    ca_coords = cmd.get_coords(f"{protein_sel} and name CA and polymer")
    if ca_coords is None or len(ca_coords) < 3:
        sys.exit(f"{protein_sel!r}: not enough CA atoms to compute an orientation")

    centroid = ca_coords.mean(axis=0)
    centered = ca_coords - centroid
    cov = centered.T @ centered
    eigvals, eigvecs = np.linalg.eigh(cov)
    long_axis = eigvecs[:, np.argmax(eigvals)]

    target = np.array([0.0, 1.0, 0.0])
    rot = rotation_to_align(long_axis, target)

    face_coords = cmd.get_coords(face_sel)
    if face_coords is not None and len(face_coords):
        face_centroid = face_coords.mean(axis=0)
        rotated_face = rot @ (face_centroid - centroid)

        if rotated_face[1] < 0:
            flip = rotation_to_align(target, -target)
            rot = flip @ rot
            rotated_face = flip @ rotated_face

        x, _, z = rotated_face
        if x * x + z * z > 1e-10:
            theta = np.arctan2(-x, z)
            c, s = np.cos(theta), np.sin(theta)
            twist = np.array([
                [c, 0, s],
                [0, 1, 0],
                [-s, 0, c],
            ])
            rot = twist @ rot

    m4 = np.eye(4)
    m4[:3, :3] = rot
    m4[:3, 3] = centroid - rot @ centroid
    cmd.transform_selection(whole_sel, list(m4.flatten()), homogenous=1)


def orient_look_down_on_plane(plane_sel, behind_sel, whole_sel):
    """Rotate `whole_sel` rigidly so the camera looks straight down onto the
    best-fit plane through `plane_sel`'s atoms (the motif hot-spot atoms,
    i.e. the protein interface) — its normal faces the camera (+Z) — rather
    than merely twisting some other axis toward it (orient_scene_vertically's
    convention, which keeps the *protein's* long axis vertical and is a poor
    fit here: the interface plane's own orientation, not the protein's
    overall long axis, is what should be face-on for a close-up of it).

    The plane normal is the smallest-variance direction of `plane_sel`'s
    atoms (PCA); its sign is chosen so `behind_sel` (the bulk of the
    protein) ends up behind the plane (-Z), not in front of it, so the
    camera looks at the interface from outside the fold. The plane's own
    largest-variance in-plane direction is made vertical (Y), giving a
    deterministic "up" without reference to the rest of the protein."""
    coords = cmd.get_coords(plane_sel)
    if coords is None or len(coords) < 3:
        sys.exit(f"{plane_sel!r}: not enough atoms to fit a plane")

    centroid = coords.mean(axis=0)
    centered = coords - centroid
    cov = centered.T @ centered
    eigvals, eigvecs = np.linalg.eigh(cov)  # ascending eigenvalue order
    normal = eigvecs[:, 0]              # smallest-variance: the plane normal
    in_plane_long_axis = eigvecs[:, 2]  # largest-variance: in-plane "long" direction

    behind_coords = cmd.get_coords(behind_sel)
    behind_centroid = behind_coords.mean(axis=0) if behind_coords is not None and len(behind_coords) else centroid
    if np.dot(behind_centroid - centroid, normal) > 0:
        normal = -normal

    v_z = normal / np.linalg.norm(normal)
    v_y = in_plane_long_axis - np.dot(in_plane_long_axis, v_z) * v_z
    v_y /= np.linalg.norm(v_y)
    v_x = np.cross(v_y, v_z)
    rot = np.array([v_x, v_y, v_z])

    m4 = np.eye(4)
    m4[:3, :3] = rot
    m4[:3, 3] = centroid - rot @ centroid
    cmd.transform_selection(whole_sel, list(m4.flatten()), homogenous=1)


def apply_material_aoshiny():
    cmd.set("ray_shadows", 0)
    cmd.set("ambient_occlusion_mode", 1)
    cmd.set("ambient_occlusion_scale", 15)
    cmd.set("ambient_occlusion_smooth", 15)
    cmd.set("ambient", 0.35)
    cmd.set("direct", 0.55)
    cmd.set("specular", 0.6)
    cmd.set("shininess", 60)
    cmd.set("spec_power", 200)
    cmd.set("spec_reflect", 1.5)
    cmd.set("reflect", 0.25)
    cmd.set("depth_cue", 0)


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


def build_hotspot_selection(protein_atoms, contig_map):
    """PyMOL selection (on the 'design' object, post-alignment) for the
    hot-spot atoms named in select_fixed_atoms, at their chai sequence
    positions. `select_fixed_atoms` names the atoms RFDiffusion3 held fixed
    during backbone generation, but the downstream LigandMPNN sequence
    design step is free to redesign a position's residue identity unless
    that position was *also* pinned for sequence design — so a motif
    residue named only by its sidechain atoms (e.g. an Arg's "NE,CD,CZ,CG,
    NH1,NH2") can still end up a different residue in the folded design.
    When none of the named atoms survive at a residue, fall back to its
    backbone atoms (always present) so the position still reads as a hot
    spot rather than silently vanishing."""
    terms = []
    for (chain, resnum), atom_names in protein_atoms.items():
        chai_pos = contig_map.get((chain, resnum))
        if chai_pos is None or not atom_names:
            continue
        present = [a for a in atom_names if cmd.count_atoms(f"design and resi {chai_pos} and name {a}") == 1]
        terms.append(f"(design and resi {chai_pos} and name {'+'.join(present or BACKBONE_ATOMS)})")
    return " or ".join(terms) if terms else "none"


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


def build_figure(panel_spec_path, out_png, protein_color, hotspot_color, ligand_color,
                  sphere_scale, ligand_representation, ligand_transparency, orient_toward,
                  zoom_buffer, trim, trim_pad, width, height, dpi, bg):
    with open(panel_spec_path) as f:
        panel = json.load(f)

    model_name = panel["model_name"]
    entry = load_design_spec(panel["design_json"], model_name)
    reference_pdb = panel.get("reference_pdb") or entry["input"]
    design_cif = panel["design_cif"]

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

    hotspot_sel = build_hotspot_selection(protein_atoms, contig_map)
    if cmd.count_atoms(hotspot_sel) == 0:
        sys.exit(f"{model_name}: hot-spot selection {hotspot_sel!r} matched no atoms in the aligned design")

    # `reference` is an M-CSA-style minimal active-site extract (the fixed
    # motif residues + ligand, ~a dozen atoms) rather than a full chain, so
    # it isn't cartoon-renderable; `design` is the actual full-length folded
    # protein, now Kabsch-aligned into the reference's coordinate frame.
    protein_sel = "design and polymer.protein"
    ligand_sel = ligand_resn_selection("reference", ligand_resns)
    if cmd.count_atoms(ligand_sel) == 0:
        sys.exit(f"{model_name}: ligand selection {ligand_sel!r} matched no atoms in {reference_pdb}")

    cmd.dss("design")
    if orient_toward == "hotspot":
        orient_look_down_on_plane(hotspot_sel, protein_sel, "reference or design")
    else:
        orient_scene_vertically(protein_sel, ligand_sel, "reference or design")

    cmd.hide("everything", "reference or design")

    cmd.show("cartoon", protein_sel)
    cmd.color(protein_color, protein_sel)

    show_ligand(ligand_sel, ligand_representation, ligand_color, ligand_transparency)

    cmd.show("spheres", hotspot_sel)
    cmd.set("sphere_scale", sphere_scale, hotspot_sel)
    cmd.color(hotspot_color, hotspot_sel)

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

    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.ray(width, height)
    cmd.png(out_png, dpi=dpi)

    if trim:
        trimmed = trim_to_content(Image.open(out_png), bg=bg, pad=trim_pad)
        trimmed.save(out_png, dpi=(dpi, dpi))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("panel_spec", help="Path to a panel-spec JSON (model_name, reference_pdb, design_json, design_cif, ...)")
    parser.add_argument("output_png")
    parser.add_argument("--protein-color", default="yellow", help="PyMOL color for the reference protein cartoon (default: yellow)")
    parser.add_argument("--hotspot-color", default="purple", help="PyMOL color for the design's motif hot-spot atom spheres (default: purple)")
    parser.add_argument("--ligand-color", default="blue", help="PyMOL color for the reference ligand (default: blue)")
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
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=1800)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--bg", default="white")
    args = parser.parse_args()

    build_figure(
        args.panel_spec, args.output_png, args.protein_color, args.hotspot_color,
        args.ligand_color, args.sphere_scale, args.ligand_representation, args.ligand_transparency,
        args.orient_toward, args.zoom_buffer, not args.no_trim, args.trim_pad,
        args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
