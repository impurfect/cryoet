"""Movie of every downloaded movie file, frame-averaged, in acquisition order.

These are the raw tilt images in cryoet_data/frames - 205 of them, 41 per tilt
series. Each .tif holds a few frames of the same view, averaged here into one.
"""
import numpy as np
import tifffile

from config import DATA, OUT
from video import binned, even, label, norm8, rgb, write

FPS = 10        # images per second
BIN = 8         # downsample factor: 5760x4092 -> 720x511

files = sorted((DATA / "frames").glob("*.tif"))
print(f"start: {len(files)} movies")

frames = []
for i, f in enumerate(files):
    img = tifffile.imread(f)
    if img.ndim == 3:
        img = img.mean(0)
    frame = rgb(norm8(binned(img, BIN)))
    frames.append(even(label(frame, f.name[:52])))
    if i % 25 == 0:
        print(f"  {i}/{len(files)}")

write(frames, OUT / "videos" / "raw_frames.mp4", FPS)
print("done")
