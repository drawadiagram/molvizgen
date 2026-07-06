#!/usr/bin/env python3
"""
GENERATE step: render the "basic problem" figure for a small-molecule binder
design campaign — just the target ligand, alone, as licorice, split into its
hot-spot atoms (the ones an RFDiffusion3 binder-design spec designates as
buried/targeted, red) and the rest of the molecule (yellow). No protein is
shown: this figure illustrates the design target, not a design.

The atom split comes from an ImpressBasePipeline `*_binder_design.json` spec
(see lib/ligand_select.py) rather than a chain id — chain-based selection
(pdz_figure.py / generate_figure.py's convention) doesn't apply here because
both atom sets live in the same ligand residue.

Usage:
    python3 ligand_hotspot_figure.py <design.json> <output.png> \
        [--task partial] [--pdb-override /path/to/ligand.pdb] \
        [--hotspot-color red] [--rest-color yellow] \
        [--width 1800] [--height 1800] [--dpi 300] [--bg white]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from ligand_select import atom_name_selection, load_binder_design_spec  # noqa: E402

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


def orient_ligand_vertically(obj_name, ligand_selection):
    """Rotate the ligand rigidly so its long axis (largest PCA component of
    all its heavy-atom coordinates) lies vertical (Y) — the single-molecule
    analog of pdz_figure.py's domain-orientation convention, for consistent
    framing across a montage even though there's no second chain to twist
    toward the camera."""
    coords = cmd.get_coords(ligand_selection)
    if coords is None or len(coords) < 3:
        sys.exit(f"{obj_name}: not enough atoms in {ligand_selection!r} to compute an orientation")

    centroid = coords.mean(axis=0)
    centered = coords - centroid
    cov = centered.T @ centered
    eigvals, eigvecs = np.linalg.eigh(cov)
    long_axis = eigvecs[:, np.argmax(eigvals)]

    target = np.array([0.0, 1.0, 0.0])
    rot = rotation_to_align(long_axis, target)

    m4 = np.eye(4)
    m4[:3, :3] = rot
    m4[:3, 3] = centroid - rot @ centroid
    cmd.transform_selection(obj_name, list(m4.flatten()), homogenous=1)


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


def build_figure(design_json, out_png, task, pdb_override, hotspot_color, rest_color,
                  width, height, dpi, bg):
    input_pdb, ligand, buried_atoms, exposed_atoms = load_binder_design_spec(design_json, task)
    pdb_path = pdb_override or input_pdb
    if not os.path.exists(pdb_path):
        sys.exit(
            f"Input structure {pdb_path!r} (from {design_json}'s {task!r} task) not found; "
            f"pass --pdb-override to point at the actual file."
        )

    obj = "ligand_target"
    cmd.load(pdb_path, obj)
    cmd.remove(f"{obj} and hydro")

    ligand_sel = f"{obj} and resn {ligand}"
    if cmd.count_atoms(ligand_sel) == 0:
        sys.exit(f"{pdb_path}: no atoms matching resn {ligand!r}")

    # Everything not in the target ligand is the "basic problem" backdrop
    # this figure omits entirely.
    cmd.remove(f"{obj} and not resn {ligand}")

    orient_ligand_vertically(obj, ligand_sel)

    cmd.hide("everything", obj)

    hotspot_sel = atom_name_selection(ligand_sel, buried_atoms)
    rest_sel = atom_name_selection(ligand_sel, exposed_atoms) if exposed_atoms else \
        f"({ligand_sel}) and not ({hotspot_sel})"

    cmd.show("sticks", hotspot_sel)
    cmd.color(hotspot_color, hotspot_sel)
    cmd.show("sticks", rest_sel)
    cmd.color(rest_color, rest_sel)
    cmd.set("stick_radius", 0.3, ligand_sel)
    cmd.set("stick_quality", 15)

    apply_material_aoshiny()

    cmd.bg_color(bg)
    cmd.zoom(obj, buffer=5)

    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.ray(width, height)
    cmd.png(out_png, dpi=dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("design_json", help="Path to an RFDiffusion3 *_binder_design.json spec")
    parser.add_argument("output_png")
    parser.add_argument("--task", default="partial", choices=["buried", "partial"],
                         help="Which task's atom split to render (default: partial, the one with a select_exposed/select_buried split)")
    parser.add_argument("--pdb-override", default=None,
                         help="Use this structure file instead of the spec's 'input' path (for when the referenced file has moved)")
    parser.add_argument("--hotspot-color", default="red", help="PyMOL color for the targeted/buried atoms (default: red)")
    parser.add_argument("--rest-color", default="yellow", help="PyMOL color for the rest of the ligand (default: yellow)")
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=1800)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--bg", default="white")
    args = parser.parse_args()

    build_figure(
        args.design_json, args.output_png, args.task, args.pdb_override,
        args.hotspot_color, args.rest_color, args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
