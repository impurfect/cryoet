"""Part 2: compare the two particle pickers.

Two pickers x two alignment branches = FOUR pick sets, and every metric below
reports all four:

  warp/etomo      warp/aretomo      pytom/etomo      pytom/aretomo

The pickers are compared WITHIN each branch, because that is the controlled
experiment: identical tomograms, only the program changes. Running it on both
branches gives two independent replications of the same comparison.

Metrics, as the assessment asks: number of detected particles, spatial overlap
between picks, detection score distributions, runtime, and key parameters.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import plotstyle
from config import (BRANCHES, DIAMETER, LABELS, MATCH_RADIUS, OUT, SIGMA,
                    SYMMETRY, TEMPLATE_EMDB, TOMO_ANGPIX)
from picks import match, pytom_picks, warp_picks
from plotstyle import COLOR, save

(OUT / "tables").mkdir(parents=True, exist_ok=True)
XYZ = ["x", "y", "z"]
SWEEP = [10, 20, 30, 40, 50, 65, 80, 100, 130, 160, 200]
VOLUME_A3 = 3470 * 4740 * 790          # tomogram volume, for the chance baseline
print("start")

warp, pytom = warp_picks(), pytom_picks()
warp["matched"] = False
pytom["matched"] = False
runtime = {"warp": json.loads((OUT / "runtime_warp_picking.json").read_text()),
           "pytom": json.loads((OUT / "runtime_pytom_picking.json").read_text())}
series = sorted(set(warp.series) & set(pytom.series))


def sphere(r):
    return 4 / 3 * np.pi * r ** 3


# ------------------------------------------- 1. counts and 2. spatial overlap
rows, pairs = [], []
for b in BRANCHES:
    for s in series:
        w = warp[(warp.branch == b) & (warp.series == s)]
        p = pytom[(pytom.branch == b) & (pytom.series == s)]
        iw, ip, dist = match(w[XYZ].to_numpy(), p[XYZ].to_numpy(), MATCH_RADIUS)
        warp.loc[w.index[iw], "matched"] = True
        pytom.loc[p.index[ip], "matched"] = True
        # what two unrelated lists of this density would agree on by accident
        chance = float(1 - np.exp(-len(p) * sphere(MATCH_RADIUS) / VOLUME_A3))
        rows.append({"branch": b, "series": s, "n_warp": len(w), "n_pytom": len(p),
                     "n_matched": len(iw),
                     "jaccard": round(len(iw) / max(len(w) + len(p) - len(iw), 1), 3),
                     "warp_confirmed": round(len(iw) / max(len(w), 1), 3),
                     "pytom_confirmed": round(len(iw) / max(len(p), 1), 3),
                     "chance": round(chance, 3),
                     "above_chance": round(len(iw) / max(len(w), 1) / max(chance, 1e-9), 2),
                     "median_separation_A": round(float(np.median(dist)), 1) if len(dist) else np.nan})
        if len(iw):
            pairs.append(pd.DataFrame({"branch": b, "series": s,
                                       "warp_score": w.iloc[iw]["score"].to_numpy(),
                                       "pytom_score": p.iloc[ip]["score"].to_numpy(),
                                       "separation": dist}))
per_tomo = pd.DataFrame(rows)
per_tomo.to_csv(OUT / "tables" / "part2_per_tomogram.csv", index=False)
pairs = pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame()

# radius sweep, per branch
sweep = []
for b in BRANCHES:
    for r in SWEEP:
        m = nw = np_ = 0
        for s in series:
            w = warp[(warp.branch == b) & (warp.series == s)][XYZ].to_numpy()
            p = pytom[(pytom.branch == b) & (pytom.series == s)][XYZ].to_numpy()
            m += len(match(w, p, r)[0]); nw += len(w); np_ += len(p)
        sweep.append({"branch": b, "radius_A": r, "radius_voxels": r / TOMO_ANGPIX,
                      "n_matched": m,
                      "jaccard": round(m / max(nw + np_ - m, 1), 3),
                      "warp_confirmed": round(m / max(nw, 1), 3)})
sweep = pd.DataFrame(sweep)
sweep.to_csv(OUT / "tables" / "part2_radius_sweep.csv", index=False)

# the same comparison with the count difference removed
eq = []
for b in BRANCHES:
    for s in series:
        w = warp[(warp.branch == b) & (warp.series == s)]
        p = pytom[(pytom.branch == b) & (pytom.series == s)]
        n = min(len(w), len(p))
        iw, _, dist = match(w.nlargest(n, "score")[XYZ].to_numpy(),
                            p.nlargest(n, "score")[XYZ].to_numpy(), MATCH_RADIUS)
        chance = float(1 - np.exp(-n * sphere(MATCH_RADIUS) / VOLUME_A3))
        eq.append({"branch": b, "series": s, "n_each": n, "n_matched": len(iw),
                   "fraction_agreeing": round(len(iw) / max(n, 1), 3),
                   "chance": round(chance, 3),
                   "above_chance": round(len(iw) / max(n, 1) / max(chance, 1e-9), 2)})
equal_n = pd.DataFrame(eq)
equal_n.to_csv(OUT / "tables" / "part2_equal_counts.csv", index=False)

# ------------------------------------------------------- 4. runtime, 4 values
rt = pd.DataFrame([{"picker": k, "branch": b, "seconds": v,
                    "seconds_per_tomogram": round(v / len(series), 1)}
                   for k, d in runtime.items() for b, v in d.items()])
rt.to_csv(OUT / "tables" / "part2_runtime.csv", index=False)

# ------------------------------------------------------------- 5. parameters
pd.DataFrame([
    {"parameter": "template", "warp": f"EMD-{TEMPLATE_EMDB}", "pytom": f"EMD-{TEMPLATE_EMDB}"},
    {"parameter": "particle diameter (A)", "warp": DIAMETER, "pytom": DIAMETER},
    {"parameter": "angular step (deg)", "warp": "7.5 (subdivisions 3)", "pytom": "7.5"},
    {"parameter": "symmetry used", "warp": f"{SYMMETRY} (octahedral, 24-fold)",
     "pytom": "C4 about z (PyTom supports z-axis symmetry only)"},
    {"parameter": "orientations searched", "warp": "1536 (36,864 / 24)",
     "pytom": "9216 (36,864 / 4)"},
    {"parameter": "search region", "warp": "auto (drops positions with <3 tilts covering them)",
     "pytom": "--search-z, set to the specimen slab"},
    {"parameter": "spectral whitening", "warp": "on", "pytom": "on"},
    {"parameter": "score definition", "warp": "sigma above background",
     "pytom": "normalised cross-correlation (LCCmax)"},
    {"parameter": "peak cutoff", "warp": f"{SIGMA} sigma",
     "pytom": f"mean + {SIGMA} x std of the correlation volume"},
]).to_csv(OUT / "tables" / "part2_parameters.csv", index=False)

# ------------------------------------------------------------------ summary
summary = []
for b in BRANCHES:
    t = per_tomo[per_tomo.branch == b]
    nw, np_, nm = int(t.n_warp.sum()), int(t.n_pytom.sum()), int(t.n_matched.sum())
    e = equal_n[equal_n.branch == b]
    pr = pairs[pairs.branch == b] if len(pairs) else pairs
    summary.append({
        "branch": b, "n_warp": nw, "n_pytom": np_, "n_matched": nm,
        "jaccard": round(nm / max(nw + np_ - nm, 1), 3),
        "warp_confirmed": round(nm / max(nw, 1), 3),
        "pytom_confirmed": round(nm / max(np_, 1), 3),
        "chance": round(float(t.chance.mean()), 3),
        "above_chance": round(float(nm / max(nw, 1) / t.chance.mean()), 2),
        "agree_at_equal_counts": round(float(e.fraction_agreeing.mean()), 3),
        "median_separation_A": round(float(pr.separation.median()), 1) if len(pr) else np.nan,
        "spearman_rho": round(float(stats.spearmanr(pr.warp_score, pr.pytom_score).statistic), 3)
        if len(pr) > 10 else np.nan,
        "warp_runtime_s": runtime["warp"][b], "pytom_runtime_s": runtime["pytom"][b]})
summary = pd.DataFrame(summary)
summary.to_csv(OUT / "tables" / "part2_summary.csv", index=False)

# are each tool's unique picks its weakest? tested per branch
uni = []
for b in BRANCHES:
    for name, df in [("warp", warp), ("pytom", pytom)]:
        d = df[df.branch == b]
        m, u = d[d.matched]["score"], d[~d.matched]["score"]
        if len(m) > 5 and len(u) > 5:
            r = stats.mannwhitneyu(m, u, alternative="greater")
            uni.append({"branch": b, "picker": name, "n_confirmed": len(m), "n_unique": len(u),
                        "median_confirmed": round(float(m.median()), 4),
                        "median_unique": round(float(u.median()), 4),
                        "prob_confirmed_higher": round(r.statistic / (len(m) * len(u)), 3),
                        "p": round(float(r.pvalue), 6)})
unique = pd.DataFrame(uni)
unique.to_csv(OUT / "tables" / "part2_unique_vs_confirmed.csv", index=False)

# =================================================================== plots
# 1. particle counts - all four sets
fig, ax = plt.subplots(figsize=(9, 3.6))
x = np.arange(len(per_tomo))
ax.bar(x - 0.27, per_tomo.n_warp, 0.25, color=COLOR["warp"], label="Warp")
ax.bar(x, per_tomo.n_matched, 0.25, color=COLOR["both"], label="agreeing")
ax.bar(x + 0.27, per_tomo.n_pytom, 0.25, color=COLOR["pytom"], label="PyTom")
ax.set_xticks(x)
ax.set_xticklabels(per_tomo.branch + "/" + per_tomo.series, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("particles detected"); ax.legend()
ax.set_title("Number of detected particles, and how many the two agree on")
save(fig, "part2_particle_counts.png")

# 2. spatial overlap - per branch, and how it depends on the tolerance
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
for b in BRANCHES:
    d = sweep[sweep.branch == b]
    axes[0].plot(d.radius_A, d.jaccard, "-o", color=COLOR[b], label=f"{LABELS[b]} tomograms")
    axes[1].plot(d.radius_A, d.warp_confirmed, "-o", color=COLOR[b], label=f"{LABELS[b]} tomograms")
for ax, ylab, title in [(axes[0], "Jaccard overlap", "Spatial overlap between pickers"),
                        (axes[1], "fraction of Warp picks confirmed",
                         "Warp picks with a PyTom pick nearby")]:
    ax.axvline(MATCH_RADIUS, color="k", ls=":", lw=1)
    ax.annotate(f"one particle\nradius ({MATCH_RADIUS} A)", (MATCH_RADIUS, 0.02),
                fontsize=6.5, xytext=(4, 0), textcoords="offset points")
    ax.set_xlabel("how close two picks must be to be the same particle (A)")
    ax.set_ylabel(ylab); ax.set_ylim(0, 1); ax.set_title(title); ax.legend(fontsize=7)
save(fig, "part2_spatial_overlap.png")

# 3. score distributions - two scales, two branches, four sets
fig, axes = plt.subplots(2, 2, figsize=(9.5, 6), sharex="col")
for row, b in enumerate(BRANCHES):
    for col, (name, df, unit) in enumerate(
            [("Warp", warp, "sigma above background"),
             ("PyTom", pytom, "normalised cross-correlation")]):
        ax = axes[row, col]
        d = df[df.branch == b]
        bins = np.linspace(df.score.min(), df.score.max(), 45)
        ax.hist(d[d.matched].score, bins=bins, alpha=0.65, color=COLOR[name.lower()],
                label=f"confirmed (n={int(d.matched.sum())})")
        ax.hist(d[~d.matched].score, bins=bins, alpha=0.65, color=COLOR["grey"],
                label=f"{name} only (n={int((~d.matched).sum())})")
        ax.legend(fontsize=6.5)
        ax.set_title(f"{name} on {LABELS[b]} tomograms", fontsize=9)
        if row == 1:
            ax.set_xlabel(f"{name} score ({unit})")
fig.suptitle("Detection score distributions - the two scores are different "
             "quantities, never a shared axis", fontsize=9)
save(fig, "part2_score_distributions.png")

# 4. runtime - all four
fig, ax = plt.subplots(figsize=(5.5, 3.4))
x = np.arange(len(BRANCHES)); wid = 0.35
for i, k in enumerate(["warp", "pytom"]):
    v = [runtime[k][b] for b in BRANCHES]
    ax.bar(x + (i - 0.5) * wid, v, wid, color=COLOR[k], label=k.capitalize())
    for xi, vi in zip(x + (i - 0.5) * wid, v):
        ax.text(xi, vi, f"{vi:.0f}s", ha="center", va="bottom", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels([f"{LABELS[b]}\ntomograms" for b in BRANCHES])
ax.set_ylabel(f"seconds, {len(series)} tomograms"); ax.legend()
ax.set_title("Picking runtime")
save(fig, "part2_runtime.png")

# 5. where the picks physically are, one tomogram per branch
fig, axes = plt.subplots(1, 2, figsize=(9, 4.6))
s0 = series[0]
for ax, b in zip(axes, BRANCHES):
    w = warp[(warp.branch == b) & (warp.series == s0)]
    p = pytom[(pytom.branch == b) & (pytom.series == s0)]
    ax.scatter(w.x, w.y, s=16, facecolors="none", edgecolors=COLOR["warp"], lw=0.7,
               label=f"Warp (n={len(w)})")
    ax.scatter(p.x, p.y, s=8, marker="x", color=COLOR["pytom"], lw=0.7,
               label=f"PyTom (n={len(p)})")
    ax.set_xlabel("x (A)"); ax.set_ylabel("y (A)")
    ax.set_aspect("equal", "datalim"); ax.legend(fontsize=7)
    ax.set_title(f"{s0} on {LABELS[b]} tomograms", fontsize=9)
fig.suptitle("Where the detected particles are, looking down the beam", fontsize=9)
save(fig, "part2_pick_positions.png")

# ---------------------------------------------------------- interpretation
lines = [
    "# Part 2 - particle-picking comparison\n",
    f"Two pickers on two alignment branches: four pick sets. The pickers are "
    f"compared within each branch, so the tomograms are identical and only the "
    f"program changes; running it on both branches replicates the comparison.\n",
    "## Summary, per branch\n", summary.to_markdown(index=False), "",
    "## Number of detected particles, per tomogram\n",
    per_tomo.to_markdown(index=False), "",
    "## Spatial overlap\n",
    f"Two picks count as the same molecule when their centres are within "
    f"{MATCH_RADIUS} A - one apoferritin radius. Beyond that the two spheres "
    f"barely overlap and calling them the same particle stops meaning anything. "
    f"Because that tolerance is the single most manipulable number here, the "
    f"full curve from {min(SWEEP)} to {max(SWEEP)} A is published alongside it.\n",
    f"`chance` is what two unrelated lists of the same density would agree on by "
    f"accident: the denser a tool's picks, the more of the other's it coincides "
    f"with for free. Without that baseline the confirmation rates are unreadable.\n",
    "## At equal counts\n", equal_n.to_markdown(index=False), "",
    f"Taking each tool's top N by score, N the smaller of the two, removes the "
    f"count difference as a confound: agreement "
    f"{equal_n.fraction_agreeing.mean():.0%} against {equal_n.chance.mean():.0%} "
    f"by chance ({equal_n.above_chance.mean():.1f}x).", "",
    "## Detection scores\n", unique.to_markdown(index=False), "",
    "`prob_confirmed_higher` is the chance a confirmed pick outscores a unique "
    "one; 0.5 means no relationship. Above 0.5 means the tool's own ranking "
    "agrees with the other tool's opinion, so its scores are informative even "
    "where its threshold is not.", "",
    "## Runtime and parameters\n", rt.to_markdown(index=False), "",
    pd.read_csv(OUT / "tables" / "part2_parameters.csv").to_markdown(index=False), "",
    "PyTom searches 9216 orientations to Warp's 1536 - a 7.5 degree search over "
    "SO(3) divided by the symmetry each can exploit. Per orientation the two run "
    "at the same speed, so the runtime gap is symmetry handling, not "
    "implementation quality.", "",
    "A higher count is not by itself a better result: a picker can return more "
    "picks purely by returning more false positives. Counts are counts, not "
    "accuracy.",
]
(OUT / "part2_interpretation.md").write_text("\n".join(lines))
print(summary.to_string(index=False))
print("done")
