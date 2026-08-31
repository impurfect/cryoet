"""Task 3: dashboard summarising both comparisons.

Everything shown here is read from results/ - no numbers or conclusions are
written into this file, so the page always reflects the current run.

Run:  streamlit run part3_dashboard.py
"""
from pathlib import Path

import pandas as pd
import streamlit as st

from config import MATCH_RADIUS, OUT, SERIES, SIGMA, TOMO_ANGPIX

st.set_page_config(page_title="Cryo-ET workflow comparison", layout="wide")
T, P, V = OUT / "tables", OUT / "plots", OUT / "videos"


def table(name):
    f = T / name
    if f.exists():
        st.dataframe(pd.read_csv(f), width="stretch", hide_index=True)
    else:
        st.info(f"{name} not found - run the analysis scripts first.")


def image(name, caption=""):
    if (P / name).exists():
        st.image(str(P / name), caption=caption, width="stretch")


def video(name, caption=""):
    if (V / name).exists():
        st.caption(caption)
        st.video(str(V / name))


st.title("Cryo-ET workflow comparison")
st.caption(f"EMPIAR-10491 apoferritin · {len(SERIES)} tilt series · "
           f"{TOMO_ANGPIX} Å/px · picks at {SIGMA}σ · match radius {MATCH_RADIUS} Å")

t1, t2, t3, t4 = st.tabs(["Alignment", "Particle picking", "Conclusions", "Data"])

with t1:
    st.header("IMOD etomo vs AreTomo2")
    st.markdown("Both branches share every step except alignment, so any "
                "difference downstream is the aligner.")
    table("part1_summary.csv")
    st.subheader("Alignment residuals")
    st.warning("Each method's own error. **Not comparable between methods** - IMOD "
               "and AreTomo minimise different quantities on different scales.")
    image("part1_alignment_residuals.png")
    st.subheader("Reconstruction quality")
    a, b = st.columns(2)
    with a:
        image("part1_reconstruction_contrast.png")
        image("part1_particles_found.png", "Downstream: particles found, and yield vs cutoff")
    with b:
        image("part1_reconstruction_sharpness.png")
        image("part1_peak_scores.png", "Downstream: how strongly molecules stood out")
    image("part1_runtime.png", "Alignment runtime")
    st.subheader("Per tilt series")
    table("part1_per_series.csv")
    video("alignment_slices.mp4", "Tomogram slices, both alignments side by side")

with t2:
    st.header("Warp vs PyTom")
    st.markdown("Both pickers ran on the same tomograms with the same template, "
                "particle diameter, angular step and peak cutoff.")
    table("part2_summary.csv")
    st.markdown("Two pickers x two alignment branches = **four pick sets**. The "
                "pickers are compared within each branch, so the tomograms are "
                "identical and only the program changes.")
    image("part2_particle_counts.png", "1. Number of detected particles")
    image("part2_spatial_overlap.png", "2. Spatial overlap, and its dependence on the tolerance")
    image("part2_score_distributions.png", "3. Detection score distributions")
    a, b = st.columns(2)
    with a:
        image("part2_runtime.png", "4. Runtime")
    with b:
        image("part2_pick_positions.png", "Where the picks are")
    st.subheader("Per tomogram")
    table("part2_per_tomogram.csv")
    st.subheader("At equal counts")
    st.caption("Each tool's top N by score, N the smaller of the two - removes the "
               "count difference as a confound.")
    table("part2_equal_counts.csv")
    st.subheader("Key parameters")
    table("part2_parameters.csv")
    st.subheader("Are the picks only one tool found its weakest?")
    table("part2_unique_vs_confirmed.csv")
    video("picks_slices.mp4", "Warp (green) and PyTom (red) picks on the same slices")

with t3:
    st.header("Conclusions")
    for f in ["part1_interpretation.md", "part2_interpretation.md"]:
        if (OUT / f).exists():
            st.markdown((OUT / f).read_text())
            st.divider()

with t4:
    st.header("The data")
    video("raw_frames.mp4", "Every raw tilt image in the dataset, in acquisition order")
    st.subheader("Radius sweep")
    table("part2_radius_sweep.csv")
