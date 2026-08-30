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
    table("task1_summary.csv")
    a, b = st.columns(2)
    with a:
        image("task1_tomogram_quality.png", "Reconstruction quality, one line per tilt series")
        image("task1_runtime.png", "Alignment runtime")
    with b:
        image("task1_particle_yield.png", "Particles found, and yield vs cutoff")
        image("task1_score_distributions.png", "Peak strength per alignment")
    st.subheader("Per tilt series")
    table("task1_per_series.csv")
    video("alignment_slices.mp4", "Tomogram slices, both alignments side by side")

with t2:
    st.header("Warp vs PyTom")
    st.markdown("Both pickers ran on the same tomograms with the same template, "
                "particle diameter, angular step and peak cutoff.")
    table("task2_summary.csv")
    a, b = st.columns(2)
    with a:
        image("task2_counts.png", "Picks per tomogram")
        image("task2_score_distributions.png", "Separate axes: the two scores are different quantities")
    with b:
        image("task2_overlap_vs_radius.png", "How agreement depends on the matching tolerance")
        image("task2_xy_example.png", "Pick positions in one tomogram")
    st.subheader("Per tomogram")
    table("task2_per_tomogram.csv")
    st.subheader("Key parameters")
    table("task2_parameters.csv")
    st.subheader("Are the picks only one tool found its weakest?")
    table("task2_unique_vs_confirmed.csv")
    video("picks_slices.mp4", "Warp (green) and PyTom (red) picks on the same slices")

with t3:
    st.header("Conclusions")
    for f in ["task1_interpretation.md", "task2_interpretation.md"]:
        if (OUT / f).exists():
            st.markdown((OUT / f).read_text())
            st.divider()

with t4:
    st.header("The data")
    video("raw_frames.mp4", "Every raw tilt image in the dataset, in acquisition order")
    st.subheader("Radius sweep")
    table("task2_radius_sweep.csv")
