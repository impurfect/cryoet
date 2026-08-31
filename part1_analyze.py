"""Part 1: compare the two alignment methods.

The assessment asks for alignment residuals, reconstruction quality measures,
runtime, and anything else informative. All four are here:

  1. alignment residuals    each aligner's own reported error (DIAGNOSTIC ONLY -
                            see the warning below)
  2. reconstruction quality contrast and sharpness of the tomograms
  3. downstream effect      particles found, and how strongly they stood out
  4. runtime

The residuals are reported because they were asked for, but they must not be
compared between methods: IMOD reports the mean distance between predicted and
observed patch positions, AreTomo reports an error from a different objective.
Different quantities, different scales. The verdict rests on 2-4, which are
measured after the branches rejoin, by the same code, on the same scale.
"""
import json
import re

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


# --------------------------------------------------------------- 1. residuals
def residuals():
    """Each aligner's own reported error, where it exposes one.

    IMOD writes a residual into its alignment logs. AreTomo does not report a
    per-series residual at all, so those rows come back empty - which is itself
    worth showing, since it is half the reason the two cannot be compared.
    """
    patterns = [
        (r"Residual error mean and sd:\s+([\d.]+)", "mean residual", "px"),
        (r"Residual error weighted mean:\s+([\d.]+)", "weighted mean residual", "px"),
        (r"Ratio of total measured values to all unknowns\s*=\s*([\d.]+)", "measured/unknowns", ""),
    ]
    cols = ["branch", "series", "metric", "value", "unit", "source"]
    rows = []
    for b, bdir in BRANCHES.items():
        dirs = sorted((bdir / "tiltstack").glob("*")) if (bdir / "tiltstack").is_dir() else []
        if not dirs:
            # AreTomo leaves no per-series log directory at all. Record that
            # rather than returning an empty frame, which would have no columns
            # and break every later reference to res.branch / res.value.
            rows.append({"branch": b, "series": "", "metric": "no logs found",
                         "value": np.nan, "unit": "", "source": ""})
            continue
        for d in dirs:
            found = False
            for log in sorted(d.glob("*.log")):
                text = log.read_text(errors="ignore")
                for pat, name, unit in patterns:
                    m = re.search(pat, text)
                    if m:
                        rows.append({"branch": b, "series": d.name, "metric": name,
                                     "value": float(m.group(1)), "unit": unit,
                                     "source": log.name})
                        found = True
            if not found:
                rows.append({"branch": b, "series": d.name, "metric": "not reported",
                             "value": np.nan, "unit": "", "source": ""})
    return pd.DataFrame(rows, columns=cols)


# ------------------------------------------------- 2. reconstruction quality
def reconstruction_quality():
    """Contrast and sharpness of each tomogram, measured at the specimen.

    The two aligners place the specimen at different depths - about 80 A apart
    on this data - so both measures locate the sample by per-slice variance
    rather than assuming it sits mid-volume.
    """
    rows = []
    for b, bdir in BRANCHES.items():
        for tomo in sorted((bdir / "reconstruction").glob("*.mrc")):
            with mrcfile.open(tomo, permissive=True) as m:
                vol = np.asarray(m.data, dtype=np.float32)
            centre = int(vol.std(axis=(1, 2)).argmax())
            half = max(vol.shape[0] // 10, 3)
            slab = vol[max(centre - half, 0):centre + half]
            rows.append({
                "branch": b, "series": tomo.stem.split("_10.00Apx")[0],
                "specimen_centre_z": centre,
                "contrast": float(slab.std() / (np.abs(slab).mean() + 1e-9)),
                "sharpness": float(ndimage.laplace(vol[centre]).var()),
            })
    return pd.DataFrame(rows)


res = residuals()
res.to_csv(OUT / "tables" / "part1_alignment_residuals.csv", index=False)
quality = reconstruction_quality()
picks = warp_picks()
runtime = json.loads((OUT / "runtime_alignment.json").read_text())

# One row per tilt series, both branches side by side.
per_series = quality.pivot(index="series", columns="branch",
                           values=["contrast", "sharpness", "specimen_centre_z"])
per_series.columns = [f"{m}_{b}" for m, b in per_series.columns]
counts = picks.groupby(["series", "branch"]).size().unstack(fill_value=0)
for b in BRANCHES:
    per_series[f"particles_{b}"] = counts[b]
    per_series[f"peak_score_{b}"] = picks[picks.branch == b].groupby("series")["score"].apply(
        lambda s: s.nlargest(50).mean())
per_series = per_series.reset_index()
per_series.to_csv(OUT / "tables" / "part1_per_series.csv", index=False)

METRICS = [("contrast", "reconstruction quality: contrast", "more structure"),
           ("sharpness", "reconstruction quality: sharpness", "crisper edges"),
           ("particles", "downstream: particles found", "more particles"),
           ("peak_score", "downstream: mean score of top 50 peaks", "stronger peaks")]

rows = []
for metric, title, higher in METRICS:
    a = per_series[f"{metric}_etomo"].to_numpy(float)
    c = per_series[f"{metric}_aretomo"].to_numpy(float)
    rows.append({"metric": title, "higher_is": higher,
                 "etomo": round(a.mean(), 5), "aretomo": round(c.mean(), 5),
                 "difference_pct": round(100 * (a - c).mean() / (abs(c.mean()) + 1e-12), 1),
                 "etomo_better_in": f"{int((a > c).sum())}/{len(a)} series",
                 "wilcoxon_p": round(float(stats.wilcoxon(a, c).pvalue), 4)})
rows.append({"metric": "alignment runtime (s, 5 series)", "higher_is": "lower is faster",
             "etomo": runtime["etomo"], "aretomo": runtime["aretomo"],
             "difference_pct": round(100 * (runtime["etomo"] - runtime["aretomo"])
                                     / runtime["aretomo"], 1),
             "etomo_better_in": "", "wilcoxon_p": ""})
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "tables" / "part1_summary.csv", index=False)

# ------------------------------------------------------------------- plots
def paired(metric, title, ylabel, name):
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for _, r in per_series.iterrows():
        ax.plot([0, 1], [r[f"{metric}_etomo"], r[f"{metric}_aretomo"]],
                "-o", color=COLOR["grey"], ms=4, lw=0.8)
        ax.annotate(r["series"], (1.03, r[f"{metric}_aretomo"]), fontsize=6, va="center")
    for i, b in enumerate(BRANCHES):
        ax.scatter([i]*len(per_series), per_series[f"{metric}_{b}"], s=50, zorder=3,
                   color=COLOR[b], label=LABELS[b])
    ax.set_xticks([0, 1]); ax.set_xticklabels(LABELS.values())
    ax.set_xlim(-0.4, 1.6); ax.set_ylabel(ylabel); ax.set_title(title)
    save(fig, name)


paired("contrast", "Reconstruction quality: contrast\n(higher = more structure)",
       "std / mean |value|", "part1_reconstruction_contrast.png")
paired("sharpness", "Reconstruction quality: sharpness\n(higher = crisper edges)",
       "variance of Laplacian", "part1_reconstruction_sharpness.png")

# residuals, per method, never on a shared axis
fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=False)
for ax, b in zip(axes, BRANCHES):
    d = res[(res.branch == b) & res.value.notna()]
    if len(d):
        ax.bar(d.series, d.value, color=COLOR[b])
        ax.set_ylabel(f"{d.metric.iloc[0]} ({d.unit.iloc[0]})")
    else:
        ax.text(0.5, 0.5, f"{LABELS[b]} reports no\nper-series residual",
                ha="center", va="center", transform=ax.transAxes, color=COLOR["grey"])
        ax.set_xticks([])
    ax.set_title(LABELS[b])
fig.suptitle("Alignment residuals - each method's own error, NOT comparable between them",
             fontsize=9)
save(fig, "part1_alignment_residuals.png")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
x = np.arange(len(per_series))
for i, b in enumerate(BRANCHES):
    axes[0].bar(x + (i - 0.5) * 0.38, per_series[f"particles_{b}"], 0.38,
                color=COLOR[b], label=LABELS[b])
axes[0].set_xticks(x); axes[0].set_xticklabels(per_series["series"])
axes[0].set_ylabel(f"particles at >= {SIGMA} sigma"); axes[0].legend()
axes[0].set_title("Downstream: particles found per tilt series")

thresholds = np.linspace(SIGMA, 12, 40)
for b in BRANCHES:
    s = picks[picks.branch == b]["score"].to_numpy()
    axes[1].plot(thresholds, [(s >= t).sum() for t in thresholds], color=COLOR[b],
                 lw=2, label=LABELS[b])
axes[1].set_xlabel("score cutoff (sigma)"); axes[1].set_ylabel("particles retained")
axes[1].set_title("Yield vs cutoff"); axes[1].legend()
save(fig, "part1_particles_found.png")

fig, ax = plt.subplots(figsize=(6, 3.4))
bins = np.linspace(SIGMA, picks["score"].max(), 50)
for b in BRANCHES:
    s = picks[picks.branch == b]["score"]
    ax.hist(s, bins=bins, alpha=0.55, color=COLOR[b], label=f"{LABELS[b]} (n={len(s)})")
ax.set_xlabel("template-matching score (sigma above background)")
ax.set_ylabel("picks"); ax.legend()
ax.set_title("Downstream: how strongly molecules stood out")
save(fig, "part1_peak_scores.png")

fig, ax = plt.subplots(figsize=(4.2, 3.4))
keys = list(BRANCHES)
ax.bar([LABELS[k] for k in keys], [runtime[k] for k in keys],
       color=[COLOR[k] for k in keys])
for i, k in enumerate(keys):
    ax.text(i, runtime[k], f"{runtime[k]:.0f}s", ha="center", va="bottom", fontsize=8)
ax.set_ylabel("seconds, 5 tilt series, 1 GPU"); ax.set_title("Alignment runtime")
save(fig, "part1_runtime.png")

# ---------------------------------------------------------- interpretation
key = summary[summary.metric.str.contains("top 50")].iloc[0]
winner = "etomo" if key["difference_pct"] > 0 else "aretomo"
faster = min(BRANCHES, key=lambda b: runtime[b])
lines = [
    "# Part 1 - alignment comparison\n", summary.to_markdown(index=False), "",
    "## Per tilt series\n", per_series.round(5).to_markdown(index=False), "",
    "## Alignment residuals\n",
    "**Reported for completeness, not for comparison.** IMOD reports the mean "
    "distance between predicted and observed tracked-patch positions. AreTomo "
    "exposes no equivalent per-series figure. Even if it did, the two minimise "
    "different objectives and are not on a common scale. Use them to spot a bad "
    "tilt series within one method, nothing more.", "",
    res.dropna(subset=["value"]).to_markdown(index=False) if res.value.notna().any()
    else "_No parseable residual found in either method's logs._", "",
    "## Verdict\n",
    f"**{LABELS[winner]}** gives the stronger downstream result: "
    f"{abs(key['difference_pct']):.1f}% higher mean score across the 50 best peaks "
    f"per tilt series ({key['etomo_better_in']} favour etomo). Reconstruction "
    f"contrast and sharpness are near-identical, so the difference shows up where "
    f"it matters - in how clearly molecules can be found.", "",
    f"**{LABELS[faster]}** is faster: {runtime[faster]:.0f} s against "
    f"{runtime['aretomo' if faster == 'etomo' else 'etomo']:.0f} s. Note this is "
    f"partly a parameter choice: AreTomo ran 5 tilt-axis refinement iterations "
    f"(--axis_iter 5) while IMOD was given the axis directly (--initial_axis).", "",
    "With five paired tilt series a Wilcoxon test cannot reach p<0.05 even when "
    "every pair agrees - 0.0625 is its floor. Consistency of direction is the "
    "evidence, not the p-value.",
]
(OUT / "part1_interpretation.md").write_text("\n".join(lines))
print(summary.to_string(index=False))
print("done")
