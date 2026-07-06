import numpy as np
import pytest

from geometry import long_axis_pca
from orient import orient_long_axis_vertical, orient_look_down_on_plane
from pymol_scene import cmd


def test_orient_long_axis_vertical_no_twist_aligns_axis_to_y(domain1_pdb):
    obj = "orient_test_no_twist"
    cmd.load(domain1_pdb, obj)
    axis_sel = f"{obj} and chain A and name CA and polymer"

    orient_long_axis_vertical(axis_sel, None, obj)

    coords = cmd.get_coords(axis_sel)
    _, long_axis = long_axis_pca(coords)
    # after orientation, the long axis should be (anti)parallel to Y
    assert abs(abs(long_axis[1]) - 1.0) < 1e-6

    cmd.delete(obj)


def test_orient_long_axis_vertical_missing_atoms_exits(domain1_pdb):
    obj = "orient_test_missing"
    cmd.load(domain1_pdb, obj)
    with pytest.raises(SystemExit):
        orient_long_axis_vertical(f"{obj} and chain Z", None, obj)
    cmd.delete(obj)


def test_orient_long_axis_vertical_twists_face_sel_above_and_toward_camera(complex_ab_pdb):
    obj = "orient_test_twist"
    cmd.load(complex_ab_pdb, obj)
    axis_sel = f"{obj} and chain A and name CA and polymer"
    face_sel = f"{obj} and chain B and name CA and polymer"

    orient_long_axis_vertical(axis_sel, face_sel, obj)

    axis_coords = cmd.get_coords(axis_sel)
    axis_centroid = axis_coords.mean(axis=0)
    _, long_axis = long_axis_pca(axis_coords)
    assert abs(abs(long_axis[1]) - 1.0) < 1e-6

    face_centroid = cmd.get_coords(face_sel).mean(axis=0)
    # face_sel ends up above the domain (+Y)...
    assert face_centroid[1] > axis_centroid[1]
    # ...and twisted to face the camera: its X/Z position should be very
    # close to the +Z axis (X near 0, Z positive).
    rel = face_centroid - axis_centroid
    assert rel[2] > 0
    assert abs(rel[0]) < 1e-6 * max(1.0, abs(rel[2])) + 1e-6

    cmd.delete(obj)


def test_orient_look_down_on_plane_puts_behind_sel_behind_camera(complex_ab_pdb):
    obj = "orient_test_plane"
    cmd.load(complex_ab_pdb, obj)
    plane_sel = f"{obj} and chain B and name CA and polymer"
    behind_sel = f"{obj} and chain A and name CA and polymer"

    orient_look_down_on_plane(plane_sel, behind_sel, obj)

    plane_centroid = cmd.get_coords(plane_sel).mean(axis=0)
    behind_centroid = cmd.get_coords(behind_sel).mean(axis=0)
    # behind_sel should end up behind the plane along -Z relative to it
    assert behind_centroid[2] < plane_centroid[2]

    cmd.delete(obj)


def test_orient_look_down_on_plane_missing_atoms_exits(domain1_pdb):
    obj = "orient_test_plane_missing"
    cmd.load(domain1_pdb, obj)
    with pytest.raises(SystemExit):
        orient_look_down_on_plane(f"{obj} and chain Z", obj, obj)
    cmd.delete(obj)
