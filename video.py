"""Shared helpers for the three video scripts."""
import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw


def read_movie(path):
    """Average the frames of one raw movie into a single image.

    These TIFFs are LZW-compressed, which tifffile can only decode with the
    imagecodecs package installed. Pillow handles LZW through libtiff and is
    already a dependency, so it stands in when imagecodecs is absent.
    """
    try:
        import tifffile
        a = tifffile.imread(path)
    except Exception:
        from PIL import ImageSequence
        a = np.stack([np.asarray(f) for f in ImageSequence.Iterator(Image.open(path))])
    return a.mean(0) if a.ndim == 3 else a


def norm8(a, lo=0.5, hi=99.5):
    """Percentile-stretch to uint8. Stretch a whole volume at once, not each
    slice, or the brightness flickers through the movie."""
    a = a.astype(np.float32)
    p0, p1 = np.percentile(a, [lo, hi])
    return (np.clip((a - p0) / (p1 - p0 + 1e-9), 0, 1) * 255).astype(np.uint8)


def binned(a, k):
    """Average k x k blocks - raw movies are far too big to show at full size."""
    h, w = a.shape[0] // k * k, a.shape[1] // k * k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))


def even(a):
    """Crop to even width and height, which h264 requires."""
    return a[:a.shape[0] // 2 * 2, :a.shape[1] // 2 * 2]


def rgb(a):
    return np.dstack([a] * 3)


def label(img, text):
    """Stamp a caption in the top-left corner of a uint8 RGB frame."""
    im = Image.fromarray(img)
    d = ImageDraw.Draw(im)
    d.rectangle([0, 0, 7 * len(text) + 8, 15], fill=(0, 0, 0))
    d.text((4, 3), text, fill=(255, 255, 255))
    return np.asarray(im)


def side_by_side(left, right):
    gap = np.zeros((left.shape[0], 4, 3), np.uint8)
    return even(np.hstack([left, gap, right]))


def write(frames, path, fps):
    path.parent.mkdir(parents=True, exist_ok=True)
    imageio.mimwrite(path, frames, fps=fps, codec="libx264", quality=8,
                     macro_block_size=1)
    print(f"wrote {path.name}  ({len(frames)} frames, {fps} fps)")
