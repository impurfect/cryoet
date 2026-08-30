"""Side-by-side movie of the two alignments' tomograms, slice by slice.

Same tilt series, same reconstruction settings - only the aligner differs, so
any difference in sharpness on screen is the alignment.
"""
import mrcfile
import numpy as np

from config import BRANCHES, LABELS, OUT, SERIES
from video import label, norm8, rgb, side_by_side, write

FPS = 10

frames = []
for s in SERIES:
    vols = {}
    for b, bdir in BRANCHES.items():
        hits = sorted((bdir / "reconstruction").glob(f"{s}_*.mrc"))
        with mrcfile.open(hits[0], permissive=True) as m:
            vols[b] = norm8(np.asarray(m.data, dtype=np.float32))
    nz = min(v.shape[0] for v in vols.values())
    for z in range(nz):
        panes = [label(rgb(vols[b][z]), f"{LABELS[b]}   {s}   z={z}") for b in BRANCHES]
        frames.append(side_by_side(*panes))
    print(f"  {s}: {nz} slices")

write(frames, OUT / "videos" / "alignment_slices.mp4", FPS)
print("done")
