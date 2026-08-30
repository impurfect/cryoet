"""Task 2: PyTom template matching, on the same tomograms Warp used.

The template comes from the same EMD-15854 map, at the voxel size in its own
header. The extraction cutoff is set to SIGMA standard deviations above each
correlation volume's background - the same rule Warp's threshold_picks applies -
so the two pickers' counts mean the same thing.
"""
import json
import subprocess
import time

import mrcfile
import numpy as np

from config import DATA, OUT, BRANCHES, TOMO_ANGPIX, DIAMETER, SIGMA, GPUS

PY = DATA / "pytom_picks"
PY.mkdir(exist_ok=True)
template, mask = PY / "template.mrc", PY / "mask.mrc"
runtime = {}

print("start")

with mrcfile.open(DATA / "emd_15854.map", permissive=True, header_only=True) as m:
    in_angpix = float(m.voxel_size.x)

subprocess.run(["pytom_create_template.py", "-i", str(DATA / "emd_15854.map"),
                "-o", str(template), "--input-voxel-size-angstrom", str(round(in_angpix, 4)),
                "--output-voxel-size-angstrom", str(TOMO_ANGPIX), "--center"], check=True)

with mrcfile.open(template, permissive=True) as m:
    box = int(m.data.shape[0])
subprocess.run(["pytom_create_mask.py", "-b", str(box), "-o", str(mask),
                "--voxel-size", str(TOMO_ANGPIX),
                "-r", str(DIAMETER / 2 / TOMO_ANGPIX), "-s", "1"], check=True)
print(f"template {box}^3 voxels from a {in_angpix:.3f} A/voxel map")

for b, bdir in BRANCHES.items():
    out = PY / b
    out.mkdir(exist_ok=True)
    t = time.time()
    for tomo in sorted((bdir / "reconstruction").glob("*.mrc")):
        name = tomo.stem.split("_10.00Apx")[0]
        subprocess.run(["pytom_match_template.py", "-t", str(template), "-m", str(mask),
                        "-v", str(tomo), "-d", str(out),
                        "--particle-diameter", str(DIAMETER), "--angular-search", "7.5",
                        "--z-axis-rotational-symmetry", "4", "--low-pass", "20",
                        "--spectral-whitening", "--random-phase-correction",
                        "-g", GPUS, "--warp-xml-file", str(bdir / f"{name}.xml")], check=True)

        with mrcfile.open(out / f"{tomo.stem}_scores.mrc", permissive=True) as m:
            d = np.asarray(m.data, dtype=np.float32)
        cutoff = float(d.mean() + SIGMA * d.std())

        subprocess.run(["pytom_extract_candidates.py",
                        "-j", str(out / f"{tomo.stem}_job.json"), "-n", "3000",
                        "--particle-diameter", str(DIAMETER),
                        "--cut-off", str(round(cutoff, 5))], check=True)
        print(f"{b}/{name}  cutoff {cutoff:.4f}")
    runtime[b] = round(time.time() - t, 1)
    print(f"{b}  {runtime[b]} s")

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "runtime_pytom_picking.json").write_text(json.dumps(runtime, indent=2))
print("done")
