"""Task 2: compare the two particle pickers.

Metrics, per the assessment: number of particles, spatial overlap between picks,
score distributions, runtime and key parameters. Both pickers ran on the same
tomograms, so the only variable is the program.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import plotstyle
from config import BRANCHES, DIAMETER, MATCH_RADIUS, OUT, SIGMA, SYMMETRY, TOMO_ANGPIX
from picks import match, pytom_picks, warp_picks
from plotstyle import COLOR, save

(OUT / "tables").mkdir(parents=True, exist_ok=True)
print("start")

warp = warp_picks()
pytom = pytom_picks()
runtime = {"warp": json.loads((OUT / "runtime_warp_picking.json").read_text()),
           "pytom": json.loads((OUT / "runtime_pytom_picking.json").read_text())}

XYZ = ["x", "y", "z"]
warp["matched"] = False
pytom["matched"] = False

# ---- counts and spatial overlap, per branch and tilt series ----
rows, pairs = [], []
for b in BRANCHES:
    for s in sorted(set(warp.series) & set(pytom.series)):
        w = warp[(warp.branch == b) & (warp.series == s)]
        p = pytom[(pytom.branch == b) & (pytom.series == s)]
        iw, ip, dist = match(w[XYZ].to_numpy(), p[XYZ].to_numpy(), MATCH_RADIUS)
        warp.loc[w.index[iw], "matched"] = True
        pytom.loc[p.index[ip], "matched"] = True
        rows.append({"branch": b, "series": s, "n_warp": len(w), "n_pytom": len(p),
                     "n_matched": len(iw),
                     "jaccard": round(len(iw) / max(len(w) + len(p) - len(iw), 1), 3),
                     "warp_confirmed": round(len(iw) / max(len(w), 1), 3),
                     "pytom_confirmed": round(len(iw) / max(len(p), 1), 3),
                     "median_separation_A": round(float(np.median(dist)), 1) if len(dist) else np.nan})
        if len(iw):
            pairs.append(pd.DataFrame({"branch": b,
                                       "warp_score": w.iloc[iw]["score"].to_numpy(),
                                       "pytom_score": p.iloc[ip]["score"].to_numpy(),
                                       "separation": dist}))
per_tomo = pd.DataFrame(rows)
per_tomo.to_csv(OUT / "tables" / "task2_per_tomogram.csv", index=False)
pairs = pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame()

# ---- how much the answer depends on the matching tolerance ----
sweep = []
for r in [10, 20, 30, 40, 50, 65, 80, 100, 130, 160, 200]:
    m = 0
    for b in BRANCHES:
        for s in sorted(set(warp.series) & set(pytom.series)):
            w = warp[(warp.branch == b) & (warp.series == s)][XYZ].to_numpy()
            p = pytom[(pytom.branch == b) & (pytom.series == s)][XYZ].to_numpy()
            m += len(match(w, p, r)[0])
    sweep.append({"radius_A": r, "radius_voxels": r / TOMO_ANGPIX, "n_matched": m,
                  "jaccard": round(m / max(len(warp) + len(pytom) - m, 1), 3)})
sweep = pd.DataFrame(sweep)
sweep.to_csv(OUT / "tables" / "task2_radius_sweep.csv", index=False)

# ---- are the picks only one tool found its weakest ones? ----
unique = []
for name, df in [("warp", warp), ("pytom", pytom)]:
    m, u = df[df.matched]["score"], df[~df.matched]["score"]
    res = stats.mannwhitneyu(m, u, alternative="greater")
    unique.append({"picker": name, "n_confirmed": len(m), "n_unique": len(u),
                   "median_confirmed": round(float(m.median()), 4),
                   "median_unique": round(float(u.median()), 4),
                   "p_confirmed_higher": round(float(res.pvalue), 6),
                   "prob_confirmed_beats_unique": round(res.statistic / (len(m) * len(u)), 3)})
unique = pd.DataFrame(unique)
unique.to_csv(OUT / "tables" / "task2_unique_vs_confirmed.csv", index=False)

rho = float(stats.spearmanr(pairs.warp_score, pairs.pytom_score).statistic) if len(pairs) else np.nan
nw, np_, nm = len(warp), len(pytom), int(per_tomo.n_matched.sum())
summary = pd.DataFrame([
    {"metric": "particles found by Warp", "value": nw},
    {"metric": "particles found by PyTom", "value": np_},
    {"metric": f"agreeing within {MATCH_RADIUS} A", "value": nm},
    {"metric": "Jaccard overlap", "value": round(nm / max(nw + np_ - nm, 1), 3)},
    {"metric": "Warp picks confirmed by PyTom", "value": round(nm / max(nw, 1), 3)},
    {"metric": "PyTom picks confirmed by Warp", "value": round(nm / max(np_, 1), 3)},
    {"metric": "median separation of agreeing picks (A)",
     "value": round(float(pairs.separation.median()), 1) if len(pairs) else np.nan},
    {"metric": "Spearman rho of scores (matched picks)", "value": round(rho, 3)},
    {"metric": "Warp runtime (s, both branches)", "value": round(sum(runtime["warp"].values()), 1)},
    {"metric": "PyTom runtime (s, both branches)", "value": round(sum(runtime["pytom"].values()), 1)},
])
summary.to_csv(OUT / "tables" / "task2_summary.csv", index=False)

pd.DataFrame([
    {"parameter": "template", "warp": f"EMD-{15854}", "pytom": f"EMD-{15854}"},
    {"parameter": "particle diameter (A)", "warp": DIAMETER, "pytom": DIAMETER},
    {"parameter": "angular step (deg)", "warp": "7.5 (subdivisions 3)", "pytom": "7.5"},
    {"parameter": "symmetry used", "warp": f"{SYMMETRY} (octahedral, 24-fold)",
     "pytom": "C4 about z (PyTom supports z-axis symmetry only)"},
    {"parameter": "spectral whitening", "warp": "on", "pytom": "on"},
    {"parameter": "score definition", "warp": "sigma above background",
     "pytom": "normalised cross-correlation (LCCmax)"},
    {"parameter": "peak cutoff", "warp": f"{SIGMA} sigma",
     "pytom": f"mean + {SIGMA} x std of the correlation volume"},
]).to_csv(OUT / "tables" / "task2_parameters.csv", index=False)

# ---- plots ----
fig, ax = plt.subplots(figsize=(8, 3.4))
lab = per_tomo.branch + "/" + per_tomo.series
x = np.arange(len(per_tomo))
ax.bar(x - 0.27, per_tomo.n_warp, 0.25, color=COLOR["warp"], label="Warp")
ax.bar(x, per_tomo.n_matched, 0.25, color=COLOR["both"], label="agreeing")
ax.bar(x + 0.27, per_tomo.n_pytom, 0.25, color=COLOR["pytom"], label="PyTom")
ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("particles"); ax.legend()
ax.set_title(f"picks per tomogram (agreement within {MATCH_RADIUS} A)")
save(fig, "task2_counts.png")

fig, ax = plt.subplots(figsize=(5.5, 3.4))
ax.plot(sweep.radius_A, sweep.jaccard, "-o", color=COLOR["both"])
ax.axvline(MATCH_RADIUS, color="k", ls=":", lw=1)
ax.annotate(f"one particle radius\n({MATCH_RADIUS} A)", (MATCH_RADIUS, 0.02),
            fontsize=7, xytext=(5, 0), textcoords="offset points")
ax.set_xlabel("how close two picks must be to count as the same particle (A)")
ax.set_ylabel("Jaccard overlap"); ax.set_ylim(0, 1)
ax.set_title("sensitivity to the matching tolerance")
save(fig, "task2_overlap_vs_radius.png")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
for ax, name, df, unit in [(axes[0], "Warp", warp, "sigma above background"),
                           (axes[1], "PyTom", pytom, "normalised cross-correlation")]:
    bins = np.linspace(df.score.min(), df.score.max(), 45)
    ax.hist(df[df.matched].score, bins=bins, alpha=0.65, color=COLOR[name.lower()],
            label=f"confirmed by the other (n={int(df.matched.sum())})")
    ax.hist(df[~df.matched].score, bins=bins, alpha=0.65, color=COLOR["grey"],
            label=f"{name} only (n={int((~df.matched).sum())})")
    ax.set_xlabel(f"{name} score ({unit})"); ax.legend(fontsize=7)
    ax.set_title(f"{name}: are unique picks the weak ones?")
save(fig, "task2_score_distributions.png")

b, s = "etomo", sorted(set(warp.series))[0]
w = warp[(warp.branch == b) & (warp.series == s)]
p = pytom[(pytom.branch == b) & (pytom.series == s)]
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(w.x, w.y, s=18, facecolors="none", edgecolors=COLOR["warp"], lw=0.8, label="Warp")
ax.scatter(p.x, p.y, s=10, marker="x", color=COLOR["pytom"], lw=0.8, label="PyTom")
ax.set_xlabel("x (A)"); ax.set_ylabel("y (A)"); ax.set_aspect("equal", "datalim")
ax.legend(); ax.set_title(f"pick positions, {b}/{s}")
save(fig, "task2_xy_example.png")

# ---- interpretation ----
probs = unique.prob_confirmed_beats_unique
verdict = ("Confirmed - for both tools the picks the other missed score lower."
           if (probs > 0.65).all() else
           "Partly confirmed - unique picks score lower but the distributions overlap heavily."
           if (probs > 0.55).all() else
           "Not confirmed - the unique picks are not reliably the low-scoring ones.")
lines = [
    "# Task 2 - particle-picking comparison\n",
    summary.to_markdown(index=False), "",
    "## Per tomogram\n", per_tomo.to_markdown(index=False), "",
    "## Key parameters\n",
    pd.read_csv(OUT / "tables" / "task2_parameters.csv").to_markdown(index=False), "",
    "## Are the picks only one tool found its weakest?\n",
    unique.to_markdown(index=False), "", verdict, "",
    f"A higher count is not by itself a better result - a picker can return more "
    f"picks purely by returning more false positives. Both tools were thresholded "
    f"by the same rule ({SIGMA} sigma above each volume's background), so the "
    f"counts are comparable, but they remain counts and not accuracy.",
    "",
    f"Agreement depends on how close two picks must be to count as the same "
    f"molecule. Across {sweep.radius_A.min()}-{sweep.radius_A.max()} A the Jaccard "
    f"overlap runs {sweep.jaccard.min():.2f} to {sweep.jaccard.max():.2f}; "
    f"{MATCH_RADIUS} A, one apoferritin radius, gives "
    f"{sweep[sweep.radius_A == MATCH_RADIUS].jaccard.iloc[0]:.2f}.",
    "",
    "Raw scores are never compared directly: Warp's is in standard deviations "
    "above background, PyTom's is a correlation coefficient. Where the two are "
    "compared it is by rank, which is scale-free.",
]
(OUT / "task2_interpretation.md").write_text("\n".join(lines))
print(summary.to_string(index=False))
print("done")
