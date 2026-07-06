import argparse
import os

import pytest

from pymol_scene import add_render_flags, apply_material_aoshiny, cmd, ray_trace_and_save


def test_apply_material_aoshiny_sets_expected_pymol_settings():
    apply_material_aoshiny()
    assert cmd.get("ray_shadows") == "off"
    assert float(cmd.get("ambient_occlusion_mode")) == 1.0
    assert cmd.get("depth_cue") == "off"
    assert float(cmd.get("ambient")) == pytest.approx(0.35)
    assert float(cmd.get("shininess")) == pytest.approx(60)


def test_add_render_flags_defaults():
    parser = argparse.ArgumentParser()
    add_render_flags(parser)
    args = parser.parse_args([])
    assert args.width == 1800
    assert args.height == 2400
    assert args.dpi == 300
    assert args.bg == "white"


def test_add_render_flags_custom_defaults():
    parser = argparse.ArgumentParser()
    add_render_flags(parser, default_height=1800)
    args = parser.parse_args([])
    assert args.height == 1800


def test_add_render_flags_overridable_from_cli():
    parser = argparse.ArgumentParser()
    add_render_flags(parser)
    args = parser.parse_args(["--width", "500", "--bg", "black"])
    assert args.width == 500
    assert args.bg == "black"


def test_ray_trace_and_save_writes_a_nonempty_png(tmp_path, domain1_pdb):
    obj = "scene_test_obj"
    cmd.load(domain1_pdb, obj)
    cmd.show("spheres", obj)
    cmd.bg_color("white")
    out_png = tmp_path / "out.png"
    ray_trace_and_save(str(out_png), width=200, height=200, dpi=72)
    cmd.delete(obj)

    assert out_png.exists()
    assert out_png.stat().st_size > 0
