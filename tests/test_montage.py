import argparse

from PIL import Image

from montage import build_montage, scale_to_width


def _save_solid(path, size, color):
    Image.new("RGB", size, color).save(path)


def test_build_montage_lays_out_grid_with_padding(tmp_path):
    p1 = tmp_path / "a.png"
    p2 = tmp_path / "b.png"
    _save_solid(p1, (100, 50), "red")
    _save_solid(p2, (100, 50), "blue")

    montage = build_montage([str(p1), str(p2)], rows=1, cols=2, padding=10, bg="white", trim=False)

    # No background to trim since images are solid-colored (no content bbox
    # differs from itself), so cell size == max image size.
    expected_w = 2 * 100 + 3 * 10
    expected_h = 1 * 50 + 2 * 10
    assert montage.size == (expected_w, expected_h)


def test_build_montage_too_many_images_for_grid_exits(tmp_path):
    p1 = tmp_path / "a.png"
    _save_solid(p1, (20, 20), "red")
    with __import__("pytest").raises(SystemExit):
        build_montage([str(p1), str(p1), str(p1)], rows=1, cols=2, padding=10, bg="white", trim=False)


def test_build_montage_centers_smaller_image_in_its_cell(tmp_path):
    small = tmp_path / "small.png"
    large = tmp_path / "large.png"
    _save_solid(small, (20, 20), "red")
    _save_solid(large, (60, 60), "blue")

    montage = build_montage([str(small), str(large)], rows=1, cols=2, padding=0, bg="white", trim=False)
    # cell size is the larger image's size (60x60), two cells wide
    assert montage.size == (120, 60)
    # small (20x20) image is centered within its (60,60) cell: offset (20,20)
    # on each side, so its content spans x,y in [20,40) within the cell.
    pixel = montage.getpixel((30, 30))
    assert pixel == (255, 0, 0)


def test_scale_to_width_preserves_aspect_ratio():
    img = Image.new("RGB", (300, 150), "white")
    scaled = scale_to_width(img, target_width_in=2.0, dpi=100)
    assert scaled.width == 200
    assert scaled.height == 100
