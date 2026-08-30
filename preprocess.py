"""Motion correction, CTF estimation and tilt-series grouping. Shared by both
alignment branches, so it runs once."""
import subprocess
from config import DATA, ANGPIX, DOSE, TOMO_DIMS, PERDEVICE


def warp(*args):
    subprocess.run(["WarpTools", *map(str, args)], cwd=DATA, check=True)


print("start")

warp("create_settings", "--folder_data", "frames", "--folder_processing", "warp_frameseries",
     "--output", "warp_frameseries.settings", "--extension", "*.tif",
     "--angpix", ANGPIX, "--gain_path", "gain_ref.mrc", "--gain_flip_y", "--exposure", DOSE)
print("frame-series settings")

warp("fs_motion_and_ctf", "--settings", "warp_frameseries.settings",
     "--m_grid", "1x1x3", "--c_grid", "2x2x1", "--c_range_max", 7, "--c_defocus_max", 8,
     "--c_use_sum", "--out_averages", "--perdevice", PERDEVICE)
print("motion + CTF")

warp("ts_import", "--mdocs", "mdoc", "--frameseries", "warp_frameseries",
     "--tilt_exposure", DOSE, "--min_intensity", 0.3, "--dont_invert", "--output", "tomostar")

warp("create_settings", "--output", "warp_tiltseries.settings",
     "--folder_processing", "warp_tiltseries", "--folder_data", "tomostar",
     "--extension", "*.tomostar", "--angpix", ANGPIX, "--gain_path", "gain_ref.mrc",
     "--gain_flip_y", "--exposure", DOSE, "--tomo_dimensions", TOMO_DIMS)
print("tilt-series grouping")

print("done")
