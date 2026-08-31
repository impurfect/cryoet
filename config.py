"""All paths and parameters. Every other script imports from here."""
import os
from pathlib import Path

# Where the raw data and Warp's processing folders live (outside the repo).
DATA = Path(os.environ.get("CRYOET_DATA", "~/cryoet_data")).expanduser()
OUT = Path(__file__).parent / "results"

# Dataset: EMPIAR-10491, the Warp tilt-series tutorial data (apoferritin).
SERIES = ["TS_1", "TS_11", "TS_17", "TS_23", "TS_32"]
ANGPIX = 0.7894            # raw pixel size, A
DOSE = 2.64                # e-/A^2 per tilt
TOMO_DIMS = "4400x6000x1000"
TOMO_ANGPIX = 10           # pixel size for alignment, reconstruction, matching
TILT_AXIS = -85.6          # initial guess, degrees

# The two alignment branches compared in Task 1.
BRANCHES = {
    "etomo": DATA / "warp_tiltseries_etomo",
    "aretomo": DATA / "warp_tiltseries_aretomo",
}
LABELS = {"etomo": "IMOD etomo", "aretomo": "AreTomo2"}

# Template: EMD-15854 is mouse apoferritin, 130 A across, octahedral symmetry.
TEMPLATE_EMDB = 15854
DIAMETER = 130             # A
SYMMETRY = "O"

SIGMA = 3                  # peak cutoff, standard deviations above background

# PyTom-specific. Warp sets its own equivalents internally.
# Resolution limit for matching, in Angstroms. 20 A is Nyquist at 10 A/voxel, so
# 20 means NO filtering - and combined with spectral whitening, which boosts high
# frequencies, that amplifies the noisiest part of the spectrum. 30 keeps the
# detail a 130 A shell actually has and drops the rest.
PYTOM_LOW_PASS = 30

# Keep only sharp, well-localised correlation peaks and reject broad diffuse
# ones. This is tuning rather than a correction, so it is off by default: report
# results with it on as "PyTom after tuning", alongside the default run.
PYTOM_TOPHAT = False
MATCH_RADIUS = 65          # A; two picks this close are the same particle
GPUS = "0"
PERDEVICE = 2
