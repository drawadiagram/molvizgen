from PIL import Image

from imgtrim import trim_to_content


def _canvas_with_square(bg="white", fg=(255, 0, 0), size=(200, 150), box=(80, 60, 120, 90)):
    img = Image.new("RGB", size, bg)
    for x in range(box[0], box[2]):
        for y in range(box[1], box[3]):
            img.putpixel((x, y), fg)
    return img


def test_trim_crops_to_content_plus_pad():
    img = _canvas_with_square()
    trimmed = trim_to_content(img, bg="white", pad=10)
    # content bbox is (80, 60, 120, 90); padded by 10 and clamped to canvas
    assert trimmed.size == (120 - 80 + 20, 90 - 60 + 20)


def test_trim_pad_clamped_to_image_bounds():
    img = _canvas_with_square(size=(200, 150), box=(0, 0, 5, 5))
    trimmed = trim_to_content(img, bg="white", pad=50)
    # pad would push left/upper negative and right/lower past bounds - clamped
    assert trimmed.size[0] <= 200
    assert trimmed.size[1] <= 150


def test_trim_all_background_returns_unchanged():
    img = Image.new("RGB", (100, 100), "white")
    trimmed = trim_to_content(img, bg="white", pad=20)
    assert trimmed is img
