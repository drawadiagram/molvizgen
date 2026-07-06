#!/usr/bin/env python3
"""
Trim uniform-background borders from a rendered figure down to its content
bounding box (plus a small pad). Every GENERATE script here ray-traces onto
a fixed-size canvas and frames the subject with a zoom buffer, so a montage
built directly from those PNGs carries a lot of incidental margin — most
visibly when panels have different aspect ratios (e.g. a tall thin ligand
next to a roughly square protein cartoon). This is the shared fix, used by
every ASSEMBLE step so panels sit close to each other regardless of how
each GENERATE script happened to frame its subject.
"""
from PIL import Image, ImageChops


def trim_to_content(image, bg="white", pad=20):
    """Crop `image` to the bounding box of everything that differs from a
    flat `bg`-colored canvas, then grow that box by `pad` pixels on every
    side (clamped to the image). Returns `image` unchanged if it's
    (near-)entirely background."""
    bg_img = Image.new(image.mode, image.size, bg)
    diff = ImageChops.difference(image, bg_img)
    bbox = diff.getbbox()
    if bbox is None:
        return image
    left, upper, right, lower = bbox
    left = max(0, left - pad)
    upper = max(0, upper - pad)
    right = min(image.width, right + pad)
    lower = min(image.height, lower + pad)
    return image.crop((left, upper, right, lower))
