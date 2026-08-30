"""Task 2: Warp's built-in 3D template matching, on both alignment branches."""
import json
import subprocess
import time
from config import DATA, OUT, BRANCHES, TOMO_ANGPIX, TEMPLATE_EMDB, DIAMETER, SYMMETRY, SIGMA

S = str(DATA / "warp_tiltseries.settings")
runtime = {}

print("start")

for b in BRANCHES:
    p = f"warp_tiltseries_{b}"
    t = time.time()
    subprocess.run(["WarpTools", "ts_template_match", "--settings", S,
                    "--input_processing", p, "--tomo_angpix", str(TOMO_ANGPIX),
                    "--subdivisions", "3", "--template_emdb", str(TEMPLATE_EMDB),
                    "--template_diameter", str(DIAMETER), "--symmetry", SYMMETRY,
                    "--whiten", "--perdevice", "1"], cwd=DATA, check=True)
    subprocess.run(["WarpTools", "threshold_picks", "--settings", S,
                    "--input_processing", p, "--in_suffix", str(TEMPLATE_EMDB),
                    "--out_suffix", "clean", "--minimum", str(SIGMA)], cwd=DATA, check=True)
    runtime[b] = round(time.time() - t, 1)
    print(f"{b}  {runtime[b]} s")

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "runtime_warp_picking.json").write_text(json.dumps(runtime, indent=2))
print("done")
