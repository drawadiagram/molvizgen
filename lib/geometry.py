#!/usr/bin/env python3
"""
Pure-numpy coordinate geometry shared by every GENERATE figure script's
orientation step: rotating one vector onto another, and extracting a
best-fit long axis or plane from a point cloud via PCA. No PyMOL dependency
here at all — the PyMOL-coupled orientation routines built on top of these
(get the coordinates, apply the rotation to a scene) live in lib/orient.py.
"""
import numpy as np


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


def _pca_axes(coords):
    """(centroid, eigvecs) for `coords` (Nx3), eigenvectors ascending by
    eigenvalue (numpy's native np.linalg.eigh order: smallest-variance
    direction first, largest-variance direction last)."""
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    cov = centered.T @ centered
    _eigvals, eigvecs = np.linalg.eigh(cov)
    return centroid, eigvecs


def long_axis_pca(coords):
    """(centroid, long_axis) for `coords` (Nx3): `long_axis` is the
    largest-variance principal axis, i.e. the point cloud's long axis."""
    if coords is None or len(coords) < 3:
        raise ValueError("long_axis_pca needs >= 3 points")
    centroid, eigvecs = _pca_axes(coords)
    return centroid, eigvecs[:, -1]


def plane_pca(coords):
    """(centroid, normal, in_plane_long_axis) for `coords` (Nx3): a best-fit
    plane through the points. `normal` is the smallest-variance principal
    axis; `in_plane_long_axis` is the largest-variance axis (necessarily
    orthogonal to `normal`), a deterministic "in-plane long direction"."""
    if coords is None or len(coords) < 3:
        raise ValueError("plane_pca needs >= 3 points")
    centroid, eigvecs = _pca_axes(coords)
    return centroid, eigvecs[:, 0], eigvecs[:, 2]
