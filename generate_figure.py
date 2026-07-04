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

Usage:
    python3 generate_figure.py <input.pdb> <output.png> \
        [--domain-chain A] [--peptide-chain B] \
        [--domain-color pink] [--peptide-color-hex 32CD32] \
        [--width 1800] [--height 2400] [--dpi 300] [--bg white]
"""
import argparse
import sys

import numpy as np

import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd  # noqa: E402


def hex_to_rgb01(hex_color):
    hex_color = hex_color.lstrip("#")
    return [int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def rotation_to_align(a, b):
    """Rotation matrix (3x3) that rotates unit vector a onto unit vector b,
    via Rodrigues' rotation formula."""
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


def orient_domain_vertically(obj_name, domain_chain, peptide_chain):
    """Rotate the whole object rigidly so the domain chain's long axis
    (largest PCA component of its CA atoms) lies along the vertical (Y)
    axis, with the peptide chain ending up above the domain and twisted
    around that vertical axis to face the camera (+Z)."""
    ca_coords = cmd.get_coords(f"{obj_name} and chain {domain_chain} and name CA and polymer")
    if ca_coords is None or len(ca_coords) < 3:
        sys.exit(f"{obj_name}: not enough chain-{domain_chain} CA atoms to compute an orientation")

    centroid = ca_coords.mean(axis=0)
    centered = ca_coords - centroid
    cov = centered.T @ centered
    eigvals, eigvecs = np.linalg.eigh(cov)
    long_axis = eigvecs[:, np.argmax(eigvals)]

    target = np.array([0.0, 1.0, 0.0])
    rot = rotation_to_align(long_axis, target)

    pep_coords = cmd.get_coords(f"{obj_name} and chain {peptide_chain} and name CA and polymer")
    if pep_coords is not None and len(pep_coords):
        pep_centroid = pep_coords.mean(axis=0)
        rotated_pep = rot @ (pep_centroid - centroid)

        if rotated_pep[1] < 0:
            flip = rotation_to_align(target, -target)
            rot = flip @ rot
            rotated_pep = flip @ rotated_pep

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


def build_figure(pdb_path, out_png, domain_chain, peptide_chain, domain_color,
                  peptide_color_hex, width, height, dpi, bg):
    obj = "complex"
    cmd.load(pdb_path, obj)
    cmd.remove(f"{obj} and hydro")

    orient_domain_vertically(obj, domain_chain, peptide_chain)

    cmd.hide("everything", obj)

    cmd.show("cartoon", f"{obj} and chain {domain_chain}")
    cmd.color(domain_color, f"{obj} and chain {domain_chain}")

    cmd.show("sticks", f"{obj} and chain {peptide_chain}")
    cmd.set_color("generate_figure_peptide_color", hex_to_rgb01(peptide_color_hex))
    cmd.set("stick_radius", 0.3, f"{obj} and chain {peptide_chain}")
    cmd.set("stick_quality", 15)
    cmd.color("generate_figure_peptide_color", f"{obj} and chain {peptide_chain}")

    apply_material_aoshiny()

    cmd.bg_color(bg)
    cmd.zoom(obj, buffer=5)

    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.ray(width, height)
    cmd.png(out_png, dpi=dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdb_path")
    parser.add_argument("output_png")
    parser.add_argument("--domain-chain", default="A", help="Chain id of the domain (default: A)")
    parser.add_argument("--peptide-chain", default="B", help="Chain id of the peptide (default: B)")
    parser.add_argument("--domain-color", default="pink", help="PyMOL color name for the domain cartoon (default: pink)")
    parser.add_argument("--peptide-color-hex", default="32CD32", help="Hex color (no #) for the peptide licorice (default: 32CD32 / limegreen)")
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=2400)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--bg", default="white")
    args = parser.parse_args()

    build_figure(
        args.pdb_path, args.output_png, args.domain_chain, args.peptide_chain,
        args.domain_color, args.peptide_color_hex, args.width, args.height, args.dpi, args.bg,
    )
    print(f"Wrote {args.output_png}")


if __name__ == "__main__":
    main()
