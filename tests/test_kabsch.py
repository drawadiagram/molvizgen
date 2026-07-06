import numpy as np
import pytest

from kabsch import as_homogeneous, kabsch_fit


def _rotation_matrix(axis, angle_deg):
    axis = axis / np.linalg.norm(axis)
    theta = np.radians(angle_deg)
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0],
    ])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def test_kabsch_fit_recovers_known_rotation_and_translation():
    rng = np.random.default_rng(0)
    mobile = rng.normal(size=(10, 3))
    R_true = _rotation_matrix(np.array([0.2, 0.7, 0.3]), 37)
    t_true = np.array([1.0, -2.0, 0.5])
    target = (mobile @ R_true.T) + t_true

    R, t = kabsch_fit(mobile, target)

    np.testing.assert_allclose(R, R_true, atol=1e-6)
    np.testing.assert_allclose(t, t_true, atol=1e-6)
    fitted = mobile @ R.T + t
    np.testing.assert_allclose(fitted, target, atol=1e-6)


def test_kabsch_fit_rejects_mismatched_or_too_few_points():
    with pytest.raises(ValueError):
        kabsch_fit(np.zeros((2, 3)), np.zeros((2, 3)))
    with pytest.raises(ValueError):
        kabsch_fit(np.zeros((5, 3)), np.zeros((4, 3)))


def test_as_homogeneous_packs_rotation_and_translation():
    R = np.eye(3) * 2
    t = np.array([1.0, 2.0, 3.0])
    flat = as_homogeneous(R, t)
    assert len(flat) == 16
    m4 = np.array(flat).reshape(4, 4)
    np.testing.assert_allclose(m4[:3, :3], R)
    np.testing.assert_allclose(m4[:3, 3], t)
    np.testing.assert_allclose(m4[3], [0, 0, 0, 1])
