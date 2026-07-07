#!/usr/bin/env python3
"""
PyMOL-coupled scene orientation, built on lib/geometry.py's pure PCA/rotation
math and lib/kabsch.py's homogeneous-matrix packing. Every GENERATE figure
script needs the whole rendered complex rotated rigidly into a consistent,
comparable framing; this module holds the two conventions used across them:

  orient_long_axis_vertical: the domain/peptide convention (pdz_figure.py,
      generate_figure.py), the ligand-only convention (ligand_hotspot_figure.py,
      no twist since there's no second thing to face the camera), and the
      arbitrary-selection convention (motif_superposition_figure.py's default
      "orient toward ligand" mode) are all the same algorithm - rotate so
      `axis_sel`'s long axis is vertical, then (if `face_sel` is given) twist
      about that vertical axis so `face_sel`'s centroid faces the camera.

  orient_look_down_on_plane: a distinct algorithm (motif_superposition_figure.py's
      "orient toward hotspot" mode) for a close-up figure that should look
      straight down onto a best-fit plane through a selection, rather than
      keep some other axis vertical.
"""
import sys

import numpy as np

from geometry import plane_pca, long_axis_pca, rotation_to_align
from kabsch import as_homogeneous
from pymol_scene import cmd


def orient_long_axis_vertical(axis_sel, face_sel, whole_sel):
    """Rotate `whole_sel` rigidly so `axis_sel`'s long axis (largest PCA
    component of its atoms) lies vertical (Y). If `face_sel` is given (not
    None and matching >0 atoms), also twist about that vertical axis so
    `face_sel`'s centroid faces the camera (+Z), flipping 180 degrees first
    if needed so it ends up above (not below) `axis_sel`'s centroid."""
    axis_coords = cmd.get_coords(axis_sel)
    try:
        centroid, long_axis = long_axis_pca(axis_coords)
    except ValueError:
        sys.exit(f"{axis_sel!r}: not enough atoms to compute an orientation")

    target = np.array([0.0, 1.0, 0.0])
    rot = rotation_to_align(long_axis, target)

    face_coords = cmd.get_coords(face_sel) if face_sel else None
    if face_coords is not None and len(face_coords):
        face_centroid = face_coords.mean(axis=0)
        rotated_face = rot @ (face_centroid - centroid)

        # Decide sign: flip 180 deg about the long axis if it puts face_sel
        # below axis_sel, so composition is consistently "axis below, face
        # above".
        if rotated_face[1] < 0:
            flip = rotation_to_align(target, -target)
            rot = flip @ rot
            rotated_face = flip @ rotated_face

        # Remaining free rotation is a twist about the (now-vertical) Y
        # axis. Choose it so face_sel's X/Z position swings to face the
        # camera (+Z), maximizing its visibility instead of leaving it
        # hidden behind axis_sel at an arbitrary twist angle.
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

    t = centroid - rot @ centroid
    cmd.transform_selection(whole_sel, as_homogeneous(rot, t), homogenous=1)


def compute_plane_orientation(plane_sel, behind_sel):
    """(rot, t) rigid transform (see lib/kabsch.py's as_homogeneous) that
    rotates so the camera looks straight down onto the best-fit plane
    through `plane_sel`'s atoms - its normal faces the camera (+Z) - rather
    than merely twisting some other axis toward it (orient_long_axis_vertical's
    convention, which keeps some other selection's long axis vertical and is
    a poor fit for a close-up of the plane itself).

    The plane normal's sign is chosen so `behind_sel` (e.g. the bulk of the
    protein) ends up behind the plane (-Z), not in front of it, so the
    camera looks at the interface from outside the fold. The plane's own
    largest-variance in-plane direction is made vertical (Y), giving a
    deterministic "up" without reference to the rest of the protein.

    Returned rather than applied directly (see orient_look_down_on_plane,
    the apply-immediately wrapper most callers want) so a transform computed
    from one selection/scene can be reapplied verbatim to a *different*
    scene sharing the same starting coordinate frame - see
    anchor_progression_zoned_figure.py, which shares one camera orientation,
    computed only from an initial/root panel's own design, across both
    panels of a two-panel figure, rather than independently re-fitting a
    plane per panel."""
    coords = cmd.get_coords(plane_sel)
    try:
        centroid, normal, in_plane_long_axis = plane_pca(coords)
    except ValueError:
        sys.exit(f"{plane_sel!r}: not enough atoms to fit a plane")

    behind_coords = cmd.get_coords(behind_sel)
    behind_centroid = behind_coords.mean(axis=0) if behind_coords is not None and len(behind_coords) else centroid
    if np.dot(behind_centroid - centroid, normal) > 0:
        normal = -normal

    v_z = normal / np.linalg.norm(normal)
    v_y = in_plane_long_axis - np.dot(in_plane_long_axis, v_z) * v_z
    v_y /= np.linalg.norm(v_y)
    v_x = np.cross(v_y, v_z)
    rot = np.array([v_x, v_y, v_z])

    t = centroid - rot @ centroid
    return rot, t


def orient_look_down_on_plane(plane_sel, behind_sel, whole_sel):
    """Rotate `whole_sel` rigidly per compute_plane_orientation(plane_sel,
    behind_sel) - the usual case, where the transform is applied immediately
    to the same scene it was computed from."""
    rot, t = compute_plane_orientation(plane_sel, behind_sel)
    cmd.transform_selection(whole_sel, as_homogeneous(rot, t), homogenous=1)
