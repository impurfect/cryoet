"""Task 1: align the same tilt series two ways, then reconstruct both.

Only the alignment step differs. Everything before it is shared (preprocess.py)
and everything after it is identical, so any difference downstream is the
aligner. Warp's --output_processing keeps the two branches apart.
"""
import json
import subprocess
import time
from config import DATA, OUT, BRANCHES, TOMO_ANGPIX, TILT_AXIS, PERDEVICE

S = str(DATA / "warp_tiltseries.settings")
runtime = {}


def warp(*args):
    subprocess.run(["WarpTools", *map(str, args)], cwd=DATA, check=True)


def timed(key, *args):
    t = time.time()
    warp(*args)
    runtime[key] = round(time.time() - t, 1)


print("start")

timed("etomo", "ts_etomo_patches", "--settings", S, "--angpix", TOMO_ANGPIX,
      "--patch_size", 500, "--initial_axis", TILT_AXIS, "--perdevice", PERDEVICE,
      "--output_processing", "warp_tiltseries_etomo")
print(f"IMOD patch tracking  {runtime['etomo']} s")

timed("aretomo", "ts_aretomo", "--settings", S, "--angpix", TOMO_ANGPIX,
      "--alignz", 800, "--axis_iter", 5, "--min_fov", 0, "--perdevice", PERDEVICE,
      "--output_processing", "warp_tiltseries_aretomo")
print(f"AreTomo2             {runtime['aretomo']} s")

for b in BRANCHES:
    p = f"warp_tiltseries_{b}"
    warp("ts_defocus_hand", "--settings", S, "--input_processing", p, "--check")
    timed(f"ctf_{b}", "ts_ctf", "--settings", S, "--input_processing", p,
          "--range_high", 7, "--defocus_max", 8, "--perdevice", PERDEVICE)
    timed(f"reconstruct_{b}", "ts_reconstruct", "--settings", S, "--input_processing", p,
          "--angpix", TOMO_ANGPIX, "--perdevice", PERDEVICE)
    print(f"CTF ({b})  {runtime[f'ctf_{b}']} s"
          f"   tomograms ({b})  {runtime[f'reconstruct_{b}']} s")

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "runtime_alignment.json").write_text(json.dumps(runtime, indent=2))
print("done")
