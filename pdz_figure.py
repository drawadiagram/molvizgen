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
import sys

import numpy as np

import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd  # noqa: E402

LIME_GREEN = [0x32 / 255, 0xCD / 255, 0x32 / 255]  # CSS "limegreen" (#32CD32)

PDZ_CHAIN = "chain A"
PEPTIDE_CHAIN = "chain B"


def rotation_to_align(a, b):
    """Rotation matrix (3x3) that rotates unit vector a onto unit vector b,
    via Rodrigues' rotation formula."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s < 1e-8:
        # a and b already parallel (or anti-parallel)
        if c > 0:
            return np.eye(3)
        # 180-degree rotation about any axis orthogonal to a
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


def orient_domain_vertically(obj_name):
    """Rotate the whole object rigidly so the PDZ domain's (chain A) long
    axis (largest PCA component of its CA atoms) lies along the vertical
    (Y) axis, with the peptide (chain B) ending up above the domain and
    twisted around that vertical axis to face the camera (+Z) so the
    peptide/binding groove isn't hidden behind the opaque surface."""
    ca_coords = cmd.get_coords(f"{obj_name} and {PDZ_CHAIN} and name CA and polymer")
    if ca_coords is None or len(ca_coords) < 3:
        sys.exit(f"{obj_name}: not enough chain-A CA atoms to compute an orientation")

    centroid = ca_coords.mean(axis=0)
    centered = ca_coords - centroid
    cov = centered.T @ centered
    eigvals, eigvecs = np.linalg.eigh(cov)
    long_axis = eigvecs[:, np.argmax(eigvals)]  # largest-variance axis

    target = np.array([0.0, 1.0, 0.0])
    rot = rotation_to_align(long_axis, target)

    pep_coords = cmd.get_coords(f"{obj_name} and {PEPTIDE_CHAIN} and name CA and polymer")
    if pep_coords is not None and len(pep_coords):
        pep_centroid = pep_coords.mean(axis=0)
        rotated_pep = rot @ (pep_centroid - centroid)

        # Decide sign: flip 180 deg about the long axis if it puts the
        # peptide below the domain, so composition is consistently
        # "domain below, peptide above".
        if rotated_pep[1] < 0:
            flip = rotation_to_align(target, -target)
            rot = flip @ rot
            rotated_pep = flip @ rotated_pep

        # Remaining free rotation is a twist about the (now-vertical) Y
        # axis. Choose it so the peptide's X/Z position swings to face the
        # camera (+Z), maximizing peptide visibility instead of leaving it
        # hidden behind the domain surface at an arbitrary twist angle.
        x, _, z = rotated_pep
        if x * x + z * z > 1e-10:
            theta = np.arctan2(-x, z)
            c, s = np.cos(theta), np.sin(theta)
            twist = np.array([
                [c, 0, s],
                [0, 1, 0],
                [-s, 0, c],
            ])
            rot = twist @ rot

    # Build a 4x4 homogeneous matrix: translate(-centroid) -> rotate -> translate(+centroid)
    m4 = np.eye(4)
    m4[:3, :3] = rot
    m4[:3, 3] = centroid - rot @ centroid
    cmd.transform_selection(obj_name, list(m4.flatten()), homogenous=1)


def apply_material_aoshiny():
    """No literal PyMOL material named 'AOShiny' exists; approximate the
    look (shiny surface + ambient occlusion, no hard shadows) via settings."""
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


def build_figure(pdb_path, out_png, width, height, dpi, bg):
    obj = "complex"
    cmd.load(pdb_path, obj)
    cmd.remove(f"{obj} and hydro")

    orient_domain_vertically(obj)

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

    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.ray(width, height)
    cmd.png(out_png, dpi=dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb_path")
    parser.add_argument("output_png")
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=2400)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--bg", default="white")
    args = parser.parse_args()

    build_figure(args.pdb_path, args.output_png, args.width, args.height, args.dpi, args.bg)
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
