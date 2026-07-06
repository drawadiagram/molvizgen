#!/usr/bin/env python3
"""
PyMOL session bootstrap and rendering boilerplate shared by every GENERATE
figure script: launching a headless PyMOL session, the "AOShiny"-style
material/lighting settings, the final ray-trace-and-save step, and the
--width/--height/--dpi/--bg flags every figure script exposes identically.

Importing this module launches the headless PyMOL session as a side effect
(same idiom lib/rmsd.py already used) — every script that needs `cmd`
should import it from here rather than repeating the finish_launching call.
"""
import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd  # noqa: E402


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


def ray_trace_and_save(out_png, width, height, dpi):
    """Ray-trace the current scene and write it to `out_png`."""
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.ray(width, height)
    cmd.png(out_png, dpi=dpi)


def add_render_flags(parser, default_width=1800, default_height=2400, default_dpi=300, default_bg="white"):
    """Add the --width/--height/--dpi/--bg flags every GENERATE figure
    script exposes identically to an existing argparse.ArgumentParser."""
    parser.add_argument("--width", type=int, default=default_width)
    parser.add_argument("--height", type=int, default=default_height)
    parser.add_argument("--dpi", type=int, default=default_dpi)
    parser.add_argument("--bg", default=default_bg)
