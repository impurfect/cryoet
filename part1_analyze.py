"""Task 1: compare the two alignments.

Metrics, per the assessment: reconstruction quality, runtime, and the downstream
effect on particle picking. The two aligners' own residuals are not compared -
IMOD and AreTomo minimise different quantities, so their numbers are not on a
common scale. Everything here is measured after the branches rejoin, where the
same code produced both.
"""
import json

import matplotlib.pyplot as plt
import mrcfile
import numpy as np
import pandas as pd
from scipy import ndimage, stats

import plotstyle
from config import BRANCHES, LABELS, OUT, SIGMA
from picks import warp_picks
from plotstyle import COLOR, save

(OUT / "tables").mkdir(parents=True, exist_ok=True)
print("start")


def tomogram_quality():
    """Contrast and sharpness of each reconstructed tomogram."""
    rows = []
    for b, bdir in BRANCHES.items():
        for tomo in sorted((bdir / "reconstruction").glob("*.mrc")):
            with mrcfile.open(tomo, permissive=True) as m:
                vol = np.asarray(m.data, dtype=np.float32)
            # Locate the specimen rather than assuming it sits mid-volume. The two
            # aligners place it at different depths - about 80 A apart on this
            # data - so measuring at nz//2 would compare one branch's sample
            # against the other's empty ice.
            per_slice = vol.std(axis=(1, 2))
            centre = int(per_slice.argmax())
            half = max(vol.shape[0] // 10, 3)
            slab = vol[max(centre - half, 0):centre + half]
            rows.append({
                "branch": b, "series": tomo.stem.split("_10.00Apx")[0],
                "sample_centre_z": centre,
                # spread of voxel values: a smeared tomogram tends to uniform grey
                "contrast": float(slab.std() / (np.abs(slab).mean() + 1e-9)),
                # variance of the Laplacian responds to edges, so it measures blur
                "sharpness": float(ndimage.laplace(vol[centre]).var()),
            })
    return pd.DataFrame(rows)


quality = tomogram_quality()
picks = warp_picks()
runtime = json.loads((OUT / "runtime_alignment.json").read_text())

# One row per tilt series, both branches side by side.
per_series = quality.pivot(index="series", columns="branch",
                           values=["contrast", "sharpness", "sample_centre_z"])
per_series.columns = [f"{m}_{b}" for m, b in per_series.columns]
counts = picks.groupby(["series", "branch"]).size().unstack(fill_value=0)
per_series["n_picks_etomo"] = counts["etomo"]
per_series["n_picks_aretomo"] = counts["aretomo"]
for b in BRANCHES:
    top = picks[picks.branch == b].groupby("series")["score"].apply(
        lambda s: s.nlargest(50).mean())
    per_series[f"top50_score_{b}"] = top
per_series = per_series.reset_index()
per_series.to_csv(OUT / "tables" / "task1_per_series.csv", index=False)

# Paired comparison: both methods saw the same five tilt series.
rows = []
for metric, higher in [("contrast", "more structure"), ("sharpness", "crisper edges"),
                       ("n_picks", "more particles"), ("top50_score", "stronger peaks")]:
    a = per_series[f"{metric}_etomo"].to_numpy(float)
    c = per_series[f"{metric}_aretomo"].to_numpy(float)
    rows.append({
        "metric": metric, "higher_is": higher,
        "etomo": round(a.mean(), 5), "aretomo": round(c.mean(), 5),
        "difference_pct": round(100 * (a - c).mean() / (abs(c.mean()) + 1e-12), 1),
        "etomo_better_in": f"{int((a > c).sum())}/{len(a)} series",
        "wilcoxon_p": round(float(stats.wilcoxon(a, c).pvalue), 4),
    })
rows.append({"metric": "alignment runtime (s)", "higher_is": "lower is faster",
             "etomo": runtime["etomo"], "aretomo": runtime["aretomo"],
             "difference_pct": round(100 * (runtime["etomo"] - runtime["aretomo"])
                                     / runtime["aretomo"], 1),
             "etomo_better_in": "", "wilcoxon_p": ""})
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "tables" / "task1_summary.csv", index=False)

# ---- plots ----
fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
for ax, metric in zip(axes, ["sharpness", "contrast"]):
    for _, r in per_series.iterrows():
        ax.plot([0, 1], [r[f"{metric}_etomo"], r[f"{metric}_aretomo"]],
                "-o", color=COLOR["grey"], ms=4, lw=0.8)
    for i, b in enumerate(BRANCHES):
        ax.scatter([i] * len(per_series), per_series[f"{metric}_{b}"],
                   s=45, zorder=3, color=COLOR[b], label=LABELS[b])
    ax.set_xticks([0, 1]); ax.set_xticklabels(LABELS.values())
    ax.set_xlim(-0.4, 1.4); ax.set_title(metric)
axes[0].set_ylabel("one line per tilt series")
save(fig, "task1_tomogram_quality.png")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
x = np.arange(len(per_series))
for i, b in enumerate(BRANCHES):
    axes[0].bar(x + (i - 0.5) * 0.38, per_series[f"n_picks_{b}"], 0.38,
                color=COLOR[b], label=LABELS[b])
axes[0].set_xticks(x); axes[0].set_xticklabels(per_series["series"])
axes[0].set_ylabel(f"particles at >= {SIGMA} sigma"); axes[0].legend()
axes[0].set_title("particles found per tilt series")

thresholds = np.linspace(SIGMA, 12, 40)
for b in BRANCHES:
    s = picks[picks.branch == b]["score"].to_numpy()
    axes[1].plot(thresholds, [(s >= t).sum() for t in thresholds],
                 color=COLOR[b], lw=2, label=LABELS[b])
axes[1].set_xlabel("score cutoff (sigma)"); axes[1].set_ylabel("particles retained")
axes[1].set_title("yield vs how strict you are"); axes[1].legend()
save(fig, "task1_particle_yield.png")

fig, ax = plt.subplots(figsize=(6, 3.4))
bins = np.linspace(SIGMA, picks["score"].max(), 50)
for b in BRANCHES:
    s = picks[picks.branch == b]["score"]
    ax.hist(s, bins=bins, alpha=0.55, color=COLOR[b], label=f"{LABELS[b]} (n={len(s)})")
ax.set_xlabel("template-matching score (sigma above background)")
ax.set_ylabel("picks"); ax.legend()
ax.set_title("how strongly molecules stood out, per alignment")
save(fig, "task1_score_distributions.png")

fig, ax = plt.subplots(figsize=(4.2, 3.4))
keys = ["etomo", "aretomo"]
ax.bar([LABELS[k] for k in keys], [runtime[k] for k in keys],
       color=[COLOR[k] for k in keys])
ax.set_ylabel("seconds, 5 tilt series, 1 GPU"); ax.set_title("alignment runtime")
save(fig, "task1_runtime.png")

# ---- interpretation, written from the numbers above ----
key = summary[summary.metric == "top50_score"].iloc[0]
winner = "etomo" if key["difference_pct"] > 0 else "aretomo"
faster = "etomo" if runtime["etomo"] < runtime["aretomo"] else "aretomo"
lines = [
    "# Task 1 - alignment comparison\n",
    summary.to_markdown(index=False), "",
    f"**{LABELS[winner]}** produces the stronger downstream result: "
    f"{abs(key['difference_pct']):.1f}% higher mean score across the 50 best peaks "
    f"per tilt series, {key['etomo_better_in'] if winner == 'etomo' else ''} "
    f"in the same direction.".replace("  ", " "),
    "",
    f"**{LABELS[faster]}** is faster: {runtime[faster]:.0f} s against "
    f"{runtime['aretomo' if faster == 'etomo' else 'etomo']:.0f} s for all five series.",
    "",
    "With five paired tilt series a Wilcoxon test cannot reach p<0.05 even when "
    "every pair agrees - 0.0625 is its floor. Consistency of direction is the "
    "evidence here, not the p-value.",
    "",
    "The two aligners' own residuals are deliberately not compared. IMOD reports "
    "the mean distance between predicted and observed patch positions; AreTomo "
    "reports an error from a different projection-matching objective. They are "
    "different quantities on different scales.",
]
(OUT / "task1_interpretation.md").write_text("\n".join(lines))
print(summary.to_string(index=False))
print("done")
