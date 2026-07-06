#!/usr/bin/env python3
"""
ASSEMBLE step: compose a wide "problem" panel occupying the left fraction of
the canvas next to an evenly-spaced rows x cols grid of "design" panels
filling the rest — the layout montage_figures.py's uniform grid can't
express on its own (a single rows x cols grid can't give one cell a
different width than the rest). Built on top of montage_figures.py: the
right-hand images are tiled with its build_montage() first, then the left
image is scaled to match that block's height and placed at the requested
width fraction of the *final* canvas, and (unless --no-scale) the whole
thing is rescaled to a target width exactly as montage_figures.py's own
final-assemble step does.

Usage:
    python3 assemble_panel_layout.py --left LEFT.png \
        --right R1.png R2.png R3.png R4.png --right-rows 2 --right-cols 2 \
        --out OUT.png [--left-width-fraction 0.3333] \
        [--target-width-in 9.0] [--dpi 300] [--padding 40] [--bg white] [--no-scale]
"""
import argparse
import os
import sys

from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from imgtrim import trim_to_content  # noqa: E402
from montage import build_montage, scale_to_width  # noqa: E402


def compose(left_path, right_paths, right_rows, right_cols, left_fraction, padding, bg, trim=True, trim_pad=20):
    right_montage = build_montage(right_paths, right_rows, right_cols, padding, bg, trim=trim, trim_pad=trim_pad)

    left_im = Image.open(left_path).convert("RGB")
    if trim:
        left_im = trim_to_content(left_im, bg, trim_pad)
    # Left panel spans the full height of the right block; width follows
    # from the requested fraction of the *total* canvas width, not of the
    # right block, so `left_fraction` reads the same way the user asked for
    # it ("the left third of the canvas").
    right_w, right_h = right_montage.width, right_montage.height
    left_w = round(right_w * left_fraction / (1 - left_fraction))

    scale = min(left_w / left_im.width, right_h / left_im.height)
    fitted_w, fitted_h = round(left_im.width * scale), round(left_im.height * scale)
    left_fitted = left_im.resize((fitted_w, fitted_h), Image.LANCZOS)

    left_cell = Image.new("RGB", (left_w, right_h), bg)
    left_cell.paste(left_fitted, ((left_w - fitted_w) // 2, (right_h - fitted_h) // 2))

    total_w = left_w + padding + right_w
    total_h = right_h
    canvas = Image.new("RGB", (total_w, total_h), bg)
    canvas.paste(left_cell, (0, 0))
    canvas.paste(right_montage, (left_w + padding, 0))
    return canvas


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--left", required=True, help="The wide left-hand panel image")
    parser.add_argument("--right", nargs="+", required=True, help="Right-hand grid images, left-to-right, top-to-bottom")
    parser.add_argument("--right-rows", type=int, default=1)
    parser.add_argument("--right-cols", type=int, default=1)
    parser.add_argument("--out", required=True, help="Output composite PNG path")
    parser.add_argument("--left-width-fraction", type=float, default=1 / 3,
                         help="Fraction of the final canvas's width the left panel occupies (default: 1/3)")
    parser.add_argument("--target-width-in", type=float, default=9.0, help="Final composite width, in inches (default: 9)")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--padding", type=int, default=40, help="Pixels between/around panels, at native resolution")
    parser.add_argument("--bg", default="white")
    parser.add_argument("--no-scale", action="store_true", help="Skip scaling to --target-width-in; save at native resolution")
    parser.add_argument("--no-trim", action="store_true", help="Don't trim each source image's background margin down to its content bounding box")
    parser.add_argument("--trim-pad", type=int, default=20, help="Pixels of background left around each image's content after trimming (default: 20)")
    args = parser.parse_args()

    if len(args.right) > args.right_rows * args.right_cols:
        sys.exit(f"{len(args.right)} right-hand images given but the {args.right_rows}x{args.right_cols} grid only has {args.right_rows * args.right_cols} cells")

    composite = compose(args.left, args.right, args.right_rows, args.right_cols, args.left_width_fraction, args.padding, args.bg,
                         trim=not args.no_trim, trim_pad=args.trim_pad)

    if args.no_scale:
        composite.save(args.out, dpi=(args.dpi, args.dpi))
        print(f"Wrote {args.out} ({composite.width}x{composite.height}px, unscaled @ {args.dpi} dpi)")
    else:
        composite = scale_to_width(composite, args.target_width_in, args.dpi)
        composite.save(args.out, dpi=(args.dpi, args.dpi))
        print(f"Wrote {args.out} ({composite.width}x{composite.height}px, {args.target_width_in}in wide @ {args.dpi} dpi)")


if __name__ == "__main__":
    main()
