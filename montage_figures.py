#!/usr/bin/env python3
"""
Cobble a set of hi-res images into a single montage: an evenly-spaced
table with the given number of (rows, columns), filled left-to-right,
top-to-bottom, then scaled down to a target width (inches, at a given DPI).

Usage:
    python3 montage_figures.py IMG1 IMG2 ... --out OUT.png \
        [--rows 1] [--cols 5] [--target-width-in 6.0] [--dpi 300] [--padding 40] [--bg white]
"""
import argparse
import sys

from PIL import Image


def build_montage(image_paths, rows, cols, padding, bg):
    if len(image_paths) > rows * cols:
        sys.exit(f"{len(image_paths)} images given but the {rows}x{cols} grid only has {rows * cols} cells")

    images = [Image.open(p).convert("RGB") for p in image_paths]

    # Evenly-spaced table: every cell is the same size (the largest image's
    # dimensions), each source image centered within its cell.
    cell_w = max(im.width for im in images)
    cell_h = max(im.height for im in images)

    total_w = cols * cell_w + (cols + 1) * padding
    total_h = rows * cell_h + (rows + 1) * padding

    canvas = Image.new("RGB", (total_w, total_h), bg)

    for idx, im in enumerate(images):
        row, col = divmod(idx, cols)  # left-to-right, top-to-bottom
        cell_x = padding + col * (cell_w + padding)
        cell_y = padding + row * (cell_h + padding)
        offset_x = cell_x + (cell_w - im.width) // 2
        offset_y = cell_y + (cell_h - im.height) // 2
        canvas.paste(im, (offset_x, offset_y))

    return canvas


def scale_to_width(image, target_width_in, dpi):
    target_px_width = round(target_width_in * dpi)
    scale = target_px_width / image.width
    target_px_height = round(image.height * scale)
    return image.resize((target_px_width, target_px_height), Image.LANCZOS)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("images", nargs="+", help="Image files, placed left-to-right, top-to-bottom")
    parser.add_argument("--rows", type=int, default=1, help="Grid rows (default: 1)")
    parser.add_argument("--cols", type=int, default=5, help="Grid columns (default: 5)")
    parser.add_argument("--out", required=True, help="Output montage PNG path")
    parser.add_argument("--target-width-in", type=float, default=6.0, help="Final montage width, in inches (default: 6)")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--padding", type=int, default=40, help="Pixels between/around cells, at native resolution")
    parser.add_argument("--bg", default="white")
    parser.add_argument("--no-scale", action="store_true", help="Skip scaling to --target-width-in; save at native padded resolution (for intermediate montages in a larger assembly)")
    args = parser.parse_args()

    montage = build_montage(args.images, args.rows, args.cols, args.padding, args.bg)
    if args.no_scale:
        montage.save(args.out, dpi=(args.dpi, args.dpi))
        print(f"Wrote {args.out} ({montage.width}x{montage.height}px, unscaled @ {args.dpi} dpi)")
    else:
        montage = scale_to_width(montage, args.target_width_in, args.dpi)
        montage.save(args.out, dpi=(args.dpi, args.dpi))
        print(f"Wrote {args.out} ({montage.width}x{montage.height}px, {args.target_width_in}in wide @ {args.dpi} dpi)")


if __name__ == "__main__":
    main()
