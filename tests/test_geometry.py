import numpy as np
import pytest

from geometry import long_axis_pca, plane_pca, rotation_to_align


def test_rotation_to_align_parallel_vectors_is_identity():
    a = np.array([1.0, 0.0, 0.0])
    rot = rotation_to_align(a, a)
    np.testing.assert_allclose(rot, np.eye(3), atol=1e-12)


def test_rotation_to_align_antiparallel_vectors_flips():
    a = np.array([1.0, 0.0, 0.0])
    rot = rotation_to_align(a, -a)
    np.testing.assert_allclose(rot @ a, -a, atol=1e-10)
    # a 180-degree rotation is its own inverse
    np.testing.assert_allclose(rot @ rot, np.eye(3), atol=1e-10)


def test_rotation_to_align_orthogonal_vectors():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    rot = rotation_to_align(a, b)
    np.testing.assert_allclose(rot @ a, b, atol=1e-10)
    # rotation matrices are orthonormal
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-10)
    np.testing.assert_allclose(np.linalg.det(rot), 1.0, atol=1e-10)


def test_rotation_to_align_general_case_maps_a_onto_b():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([-2.0, 0.5, 1.0])
    rot = rotation_to_align(a, b)
    mapped = rot @ (a / np.linalg.norm(a))
    np.testing.assert_allclose(mapped, b / np.linalg.norm(b), atol=1e-10)


def test_long_axis_pca_finds_dominant_elongation_direction():
    rng = np.random.default_rng(1)
    n = 200
    # elongated along X, small jitter in Y/Z
    coords = np.column_stack([
        np.linspace(-10, 10, n),
        rng.normal(scale=0.1, size=n),
        rng.normal(scale=0.1, size=n),
    ])
    centroid, long_axis = long_axis_pca(coords)
    np.testing.assert_allclose(centroid, [0, 0, 0], atol=0.1)
    # long_axis should be (anti)parallel to the X axis
    assert abs(abs(long_axis[0]) - 1.0) < 1e-3
    assert abs(long_axis[1]) < 1e-2
    assert abs(long_axis[2]) < 1e-2


def test_long_axis_pca_requires_at_least_3_points():
    with pytest.raises(ValueError):
        long_axis_pca(np.zeros((2, 3)))
    with pytest.raises(ValueError):
        long_axis_pca(None)


def test_plane_pca_finds_normal_and_in_plane_axis():
    rng = np.random.default_rng(2)
    n = 200
    # points scattered widely in X (most variance), moderately in Y, and
    # tightly clustered in Z (least variance) -> normal should be ~Z axis.
    coords = np.column_stack([
        rng.normal(scale=5.0, size=n),
        rng.normal(scale=2.0, size=n),
        rng.normal(scale=0.05, size=n),
    ])
    centroid, normal, in_plane_long_axis = plane_pca(coords)
    np.testing.assert_allclose(centroid, [0, 0, 0], atol=0.5)
    assert abs(abs(normal[2]) - 1.0) < 1e-2
    assert abs(abs(in_plane_long_axis[0]) - 1.0) < 1e-2
    # normal and in-plane axis are orthogonal
    assert abs(np.dot(normal, in_plane_long_axis)) < 1e-8


def test_plane_pca_requires_at_least_3_points():
    with pytest.raises(ValueError):
        plane_pca(np.zeros((1, 3)))
