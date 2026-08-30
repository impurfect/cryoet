"""Load both pickers' output into one tidy table, in Angstroms.

Warp writes coordinates NORMALISED to 0-1 across the volume and scores its peaks
in standard deviations above background. PyTom writes VOXELS and a correlation
coefficient. Both are converted to Angstroms here so positions can be compared;
the scores are left on their own scales and are never mixed.
"""
import mrcfile
import numpy as np
import pandas as pd
import starfile

from config import BRANCHES, DATA, TOMO_ANGPIX


def _shape(branch_dir, series):
    hits = list((branch_dir / "reconstruction").glob(f"{series}_*.mrc"))
    with mrcfile.open(hits[0], permissive=True) as m:
        nz, ny, nx = m.data.shape
    return np.array([nx, ny, nz], dtype=float)


def _read(path):
    df = starfile.read(path)
    return next(iter(df.values())) if isinstance(df, dict) else df


def warp_picks():
    rows = []
    for b, bdir in BRANCHES.items():
        for star in sorted((bdir / "matching").glob("*clean.star")):
            series = star.stem.split("_10.00Apx")[0]
            df = _read(star)
            if len(df) == 0:
                continue
            xyz = df[["rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]].to_numpy(float)
            xyz = xyz * _shape(bdir, series) * TOMO_ANGPIX      # 0-1 -> A
            rows.append(pd.DataFrame({
                "branch": b, "series": series, "picker": "warp",
                "x": xyz[:, 0], "y": xyz[:, 1], "z": xyz[:, 2],
                "score": df["rlnAutopickFigureOfMerit"].to_numpy(float)}))
    return pd.concat(rows, ignore_index=True)


def pytom_picks():
    rows = []
    for b in BRANCHES:
        for star in sorted((DATA / "pytom_picks" / b).glob("*particles.star")):
            series = star.stem.split("_10.00Apx")[0]
            df = _read(star)
            if len(df) == 0:
                continue
            xyz = df[["rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]].to_numpy(float)
            rows.append(pd.DataFrame({
                "branch": b, "series": series, "picker": "pytom",
                "x": xyz[:, 0] * TOMO_ANGPIX, "y": xyz[:, 1] * TOMO_ANGPIX,
                "z": xyz[:, 2] * TOMO_ANGPIX,
                "score": df["rlnLCCmax"].to_numpy(float)}))
    return pd.concat(rows, ignore_index=True)


def match(a, b, radius):
    """Pair two sets of 3D positions, each point used at most once.

    Solved as an assignment problem (Hungarian algorithm) rather than by nearest
    neighbour, so no point is double counted and the answer does not depend on
    the order rows happen to appear in the file.
    """
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.distance import cdist
    if len(a) == 0 or len(b) == 0:
        return np.array([], int), np.array([], int), np.array([], float)
    d = cdist(a, b)
    rows, cols = linear_sum_assignment(np.where(d <= radius, d, radius * 1e6))
    keep = d[rows, cols] <= radius
    return rows[keep], cols[keep], d[rows, cols][keep]
