#!/usr/bin/env python3
"""
Image-tiling and resizing helpers shared by the two ASSEMBLE steps:
montage_figures.py (an evenly-spaced rows x cols grid) and
assemble_panel_layout.py (a wide left panel next to such a grid, built by
calling build_montage() directly for its right-hand block).
"""
import sys

from PIL import Image

from imgtrim import trim_to_content


def build_montage(image_paths, rows, cols, padding, bg, trim=True, trim_pad=20):
    if len(image_paths) > rows * cols:
        sys.exit(f"{len(image_paths)} images given but the {rows}x{cols} grid only has {rows * cols} cells")

    images = [Image.open(p).convert("RGB") for p in image_paths]
    if trim:
        images = [trim_to_content(im, bg, trim_pad) for im in images]

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
