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


def render_solo(obj, other_objs, out_png, zoom_buffer, width, height, dpi, bg):
    """Ray-trace `obj` alone out of a PyMOL session that also holds
    `other_objs` (e.g. a Kabsch-aligned reference/design pair sharing one
    coordinate frame): disable every name in `other_objs` so they don't
    appear in this render, zoom to fit just `obj`, ray-trace, then re-enable
    them so the session's visibility state doesn't leak into the next
    render this script call makes (a combined/overlaid view, or the other
    object's own solo view).

    Any figure script that renders a combined/overlaid view of multiple
    objects should offer this as an option alongside it (see
    aligned_overlay_figure.py's --out-reference/--out-design and
    aligned_pair_figure.py, which is built entirely out of two calls to
    this) -- a reader comparing two structures often wants each one
    unobstructed, not just the overlay."""
    for name in other_objs:
        cmd.disable(name)
    cmd.enable(obj)
    cmd.bg_color(bg)
    # cmd.zoom fits to the *current* viewport's aspect ratio, which
    # defaults to a headless session's 640x480 (4:3) regardless of the
    # width/height the figure is eventually ray-traced at -- for a portrait
    # canvas (the default 1800x2400, 3:4) that mismatch crops/off-centers
    # whatever was just zoomed. Set the real output aspect ratio first so
    # the zoom is computed for the frame it will actually be rendered into.
    cmd.viewport(width, height)
    cmd.zoom(obj, buffer=zoom_buffer)
    ray_trace_and_save(out_png, width, height, dpi)
    for name in other_objs:
        cmd.enable(name)


def add_render_flags(parser, default_width=1800, default_height=2400, default_dpi=300, default_bg="white"):
    """Add the --width/--height/--dpi/--bg flags every GENERATE figure
    script exposes identically to an existing argparse.ArgumentParser."""
    parser.add_argument("--width", type=int, default=default_width)
    parser.add_argument("--height", type=int, default=default_height)
    parser.add_argument("--dpi", type=int, default=default_dpi)
    parser.add_argument("--bg", default=default_bg)
