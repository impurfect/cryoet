"""Side-by-side movie of the same tomogram slices, annotated with each picker's
detections: Warp in green, PyTom in red.

A pick is drawn on a slice when its centre lies within one particle radius of
that slice, so each particle appears across the few slices it actually spans.
"""
import mrcfile
import numpy as np
from PIL import Image, ImageDraw

from config import BRANCHES, DIAMETER, OUT, SERIES, TOMO_ANGPIX
from picks import pytom_picks, warp_picks
from video import label, norm8, rgb, side_by_side, write

FPS = 10
BRANCH = "etomo"                       # tomograms held constant for Task 2
R = DIAMETER / 2 / TOMO_ANGPIX         # particle radius in voxels

sets = [("Warp", warp_picks(), (60, 255, 60)),
        ("PyTom", pytom_picks(), (255, 70, 70))]
print("start")

frames = []
for s in SERIES:
    hits = sorted((BRANCHES[BRANCH] / "reconstruction").glob(f"{s}_*.mrc"))
    with mrcfile.open(hits[0], permissive=True) as m:
        vol = norm8(np.asarray(m.data, dtype=np.float32))
    picks = [(n, df[(df.branch == BRANCH) & (df.series == s)], c) for n, df, c in sets]

    for z in range(vol.shape[0]):
        panes = []
        for name, df, colour in picks:
            near = df[np.abs(df.z / TOMO_ANGPIX - z) <= R]
            im = Image.fromarray(rgb(vol[z]))
            draw = ImageDraw.Draw(im)
            for x, y in zip(near.x / TOMO_ANGPIX, near.y / TOMO_ANGPIX):
                draw.ellipse([x - R, y - R, x + R, y + R], outline=colour, width=1)
            panes.append(label(np.asarray(im), f"{name}   {s}   z={z}   n={len(near)}"))
        frames.append(side_by_side(*panes))
    print(f"  {s}: {vol.shape[0]} slices")

write(frames, OUT / "videos" / "picks_slices.mp4", FPS)
print("done")
