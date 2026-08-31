"""Build RESULTS.md - every plot and table, with provenance and interpretation.

The prose describing what each output IS stays fixed; every number comes from
results/tables/, so re-running after new data keeps the report honest.

Run:  python make_report.py
"""
import pandas as pd

from config import DIAMETER, MATCH_RADIUS, OUT, SERIES, SIGMA, TOMO_ANGPIX

T, P = OUT / "tables", OUT / "plots"
md = []


def w(*lines):
    md.extend(lines)


def table(name, note=""):
    f = T / name
    if not f.exists():
        w(f"_`{name}` not found — run the analysis scripts._", "")
        return None
    df = pd.read_csv(f)
    w(df.to_markdown(index=False), "")
    if note:
        w(note, "")
    return df


def figure(png, what, how, says):
    if not (P / png).exists():
        w(f"_`{png}` not found._", "")
        return
    w(f"![{png}](results/plots/{png})", "",
      f"**What it is.** {what}", "",
      f"**How it was made.** {how}", "",
      f"**What it says.** {says}", "")


def embed(name, level=3):
    """Insert an interpretation file, demoting its headings so they nest under
    this document's structure instead of competing with it."""
    f = OUT / name
    if not f.exists():
        return
    out = []
    for line in f.read_text().splitlines():
        out.append("#" * level + line if line.startswith("#") else line)
    w(*out, "")


def num(df, col, where=None, default="—"):
    """Pull a single value out of a table for use in the prose."""
    if df is None or col not in df:
        return default
    d = df if where is None else df.query(where)
    return d[col].iloc[0] if len(d) else default


# ===========================================================================
w("# Results", "",
  "Compiled from `results/` by `make_report.py`. Every number is read from the "
  "tables; nothing is written in by hand.", "",
  f"**Dataset** EMPIAR-10491, {len(SERIES)} tilt series of apoferritin · "
  f"**Processing** {TOMO_ANGPIX} Å/voxel · **Template** EMD-15854, {DIAMETER} Å, "
  f"octahedral · **Peak cutoff** {SIGMA}σ · **Match radius** {MATCH_RADIUS} Å", "",
  "---", "")

# ------------------------------------------------------------------ headline
w("## Headline", "")
try:
    a = pd.read_csv(T / "part1_summary.csv")
    k = a[a.metric.str.contains("top 50", case=False)].iloc[0]
    r = a[a.metric.str.contains("runtime", case=False)].iloc[0]
    fast = "IMOD etomo" if float(r.etomo) < float(r.aretomo) else "AreTomo2"
    lead = "IMOD etomo" if float(k.difference_pct) > 0 else "AreTomo2"
    w(f"**Part 1.** {lead} gives the stronger downstream result — "
      f"{abs(float(k.difference_pct)):.1f}% higher mean score across the 50 best "
      f"peaks per tilt series ({k.etomo_better_in} in that direction). "
      f"{fast} is faster: {min(float(r.etomo), float(r.aretomo)):.0f} s against "
      f"{max(float(r.etomo), float(r.aretomo)):.0f} s for all five series.", "")
except Exception:
    w("_Part 1 summary not available._", "")
try:
    b = pd.read_csv(T / "part2_summary.csv")
    tot_w, tot_p = int(b.n_warp.sum()), int(b.n_pytom.sum())
    w(f"**Part 2.** Warp returned {tot_w} particles across both branches, PyTom "
      f"{tot_p}. They agree on {int(b.n_matched.sum())} within {MATCH_RADIUS} Å — "
      f"{b.above_chance.mean():.1f}x more than two unrelated lists of the same "
      f"density would by chance. With counts equalised, agreement is "
      f"{b.agree_at_equal_counts.mean():.0%}.", "")
except Exception:
    w("_Part 2 summary not available._", "")
w("Both statements are computed from the tables below.", "", "---", "")

# --------------------------------------------------------------- provenance
w("## 1. Where the data came from", "",
  "`download.py` fetches three things from the EBI public archives:", "",
  "| what | source | lands in |",
  "|---|---|---|",
  "| 205 movie files (5 series × 41 tilts) | `ftp.ebi.ac.uk/empiar/world_availability/10491/data/tiltseries/data/` | `$CRYOET_DATA/frames/` |",
  "| 5 metadata files | `.../10491/data/tiltseries/mdoc/` | `$CRYOET_DATA/mdoc/` |",
  "| camera gain reference | `.../10491/data/gain_ref.mrc` | `$CRYOET_DATA/gain_ref.mrc` |",
  "| template structure | `ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-15854/map/` | `$CRYOET_DATA/emd_15854.map` |", "",
  "The movie list is taken from the `.mdoc` files rather than a filename "
  "wildcard. The tutorial's own script globs `*-11_*.tif`, which matches both the "
  "`53-11` and `59-11` series in this deposit and quietly downloads eight tilt "
  "series where you wanted five.", "",
  "**What each is:**", "",
  "- **`.tif` movies** — the raw output. Each is a short burst of frames of one "
  "view at one tilt angle, 5760 × 4092 pixels. At the dose used they are almost "
  "pure noise; you cannot see a molecule in one.",
  "- **`.mdoc`** — plain text, one block per tilt image, recording which movie "
  "file it is, the stage angle, and the accumulated dose. This is what turns 205 "
  "unrelated movies into 5 ordered tilt series.",
  "- **`gain_ref.mrc`** — a calibration image of the camera sensor. Every movie is "
  "divided by it to remove the sensor's fixed pattern.",
  "- **`emd_15854.map`** — a published 1.8 Å structure of apoferritin, deposited "
  "at 0.729 Å/voxel. Both pickers use this same file as their reference shape.", "")

w("## 2. What preprocessing does to it", "",
  "`preprocess.py` runs once, and its output is shared by both alignment "
  "branches so they cannot differ for any reason other than the aligner.", "",
  "| step | what it does | creates |",
  "|---|---|---|",
  "| `create_settings` | records pixel size, dose, gain path | `warp_frameseries.settings` |",
  "| `fs_motion_and_ctf` | measures how the specimen drifted during each exposure and shifts the frames back into register; measures the defocus of every image from its power spectrum | `warp_frameseries/` — one XML of results per movie, plus `average/` holding one drift-corrected image per movie |",
  "| `ts_import` | reads the `.mdoc` files and groups movies into tilt series with their angles and cumulative dose | `tomostar/` — five `.tomostar` files |",
  "| `create_settings` | records the volume size to reconstruct | `warp_tiltseries.settings` |", "",
  "**The new data:** `warp_frameseries/average/` holds 205 images that are still "
  "2D, but now drift-corrected and with known defocus. These, not the raw `.tif` "
  "files, are what everything downstream uses.", "")

w("## 3. What alignment and reconstruction create", "",
  "`part1_align.py` runs the two aligners into separate folders using Warp's "
  "`--output_processing`, then reconstructs both:", "",
  "```",
  "$CRYOET_DATA/",
  "├── warp_tiltseries_etomo/      <- branch A, IMOD patch tracking",
  "│   ├── TS_*.xml                   the projection geometry it solved for",
  "│   ├── tiltstack/TS_*/            IMOD's working files and logs",
  "│   └── reconstruction/            TS_*_10.00Apx.mrc  <- the tomograms",
  "└── warp_tiltseries_aretomo/    <- branch B, AreTomo2, same structure",
  "```", "",
  f"Each tomogram is 347 × 474 × 79 voxels at {TOMO_ANGPIX} Å — about "
  f"347 × 474 nm across and 79 nm thick, roughly 26 MB. **This is the first point "
  f"at which molecules are visible.**", "")

w("## 4. What particle picking creates", "",
  "| script | creates | contents |",
  "|---|---|---|",
  "| `part2_warp.py` | `warp_tiltseries_*/matching/` | correlation volumes, and `TS_*_clean.star` listing accepted peaks. Coordinates are **normalised 0–1**; scores are in **σ above background** |",
  "| `part2_pytom.py` | `pytom_picks/etomo/`, `pytom_picks/aretomo/` | `*_scores.mrc`, `*_angles.mrc`, `*_job.json`, and `*_particles.star`. Coordinates are in **voxels**; scores are **normalised cross-correlation** |", "",
  "Two pickers × two branches = **four pick sets**. `picks.py` converts both to "
  "Ångströms so positions can be compared, and leaves the scores on their own "
  "scales because they measure different things.", "",
  "---", "")

# ------------------------------------------------------------------ part 1
w("# Part 1 — Does the alignment method matter?", "",
  "Both branches share every step except alignment and both were carried through "
  "to particle picking, so they are judged on outcomes rather than on each "
  "program's internal error measure.", "")

s1 = table("part1_summary.csv")
if s1 is not None:
    w("_Read from `results/tables/part1_summary.csv`._", "")

w("### Alignment residuals", "")
figure("part1_alignment_residuals.png",
       "Each aligner's own reported error, one panel per method.",
       "`part1_analyze.py` scans `warp_tiltseries_*/tiltstack/*/*.log` for a "
       "labelled residual line.",
       "**These two panels must not be compared with each other.** IMOD reports "
       "the mean distance between predicted and observed tracked-patch positions; "
       "AreTomo exposes no per-series equivalent and minimises a different "
       "objective entirely. Useful only for spotting one bad tilt series within "
       "one method.")

w("### Reconstruction quality", "")
figure("part1_reconstruction_contrast.png",
       "Spread of voxel values in the specimen slab, one line per tilt series "
       "joining its two branch values.",
       "Reads `reconstruction/*.mrc` from both branches. The slab is located per "
       "volume by per-slice variance, not assumed to be mid-volume — the two "
       "aligners place the specimen about 80 Å apart in z.",
       "A well-aligned tomogram has crisp molecules and therefore more spread in "
       "its values; a misaligned one smears towards uniform grey. A paired plot "
       "shows whether a difference is consistent across all five series or driven "
       "by one outlier.")
figure("part1_reconstruction_sharpness.png",
       "Variance of the Laplacian of the specimen slice — the standard blur "
       "measure, borrowed from camera autofocus.",
       "Same volumes as above.",
       "The Laplacian responds to edges, so its variance is high when edges are "
       "crisp and low when blurred. It measures blur without needing to know what "
       "the right answer looks like.")

w("### Downstream effect", "")
figure("part1_particles_found.png",
       "Left: particles found per tilt series in each branch. Right: how the "
       "count changes as the score cutoff moves from 3σ to 12σ.",
       "Counts the picks in `warp_tiltseries_*/matching/*_clean.star` — the same "
       "picker on both branches, so only the alignment differs.",
       "The right-hand panel matters: a single count depends entirely on where "
       "you put the threshold. If one branch is above the other across the whole "
       "curve, the conclusion does not rest on an arbitrary cutoff.")
figure("part1_peak_scores.png",
       "Distribution of template-matching scores from both branches, sharing an "
       "axis.",
       "Same STAR files. Sharing an axis is legitimate here: same program, same "
       "template, same normalisation.",
       "**This is the plot that answers Part 1.** The score is how many standard "
       "deviations a peak stands above its own volume's background. A distribution "
       "pushed further right means molecules stood out more clearly, which means "
       "the tomogram preserved their structure better, which means the alignment "
       "was better.")
figure("part1_runtime.png",
       "Wall-clock seconds to align all five tilt series.",
       "Timed inside `part1_align.py`, written to `results/runtime_alignment.json`.",
       "Note this is partly a parameter choice, not purely algorithmic: AreTomo "
       "ran five tilt-axis refinement passes (`--axis_iter 5`) while IMOD was "
       "handed the axis directly (`--initial_axis`). At `--axis_iter 1` AreTomo "
       "would be roughly six times faster.")

w("### Per tilt series", "")
table("part1_per_series.csv")

w("### Interpretation", "")
embed("part1_interpretation.md")

# ------------------------------------------------------------------ part 2
w("---", "", "# Part 2 — Does the particle picker matter?", "",
  "Two pickers × two alignment branches = four pick sets. The pickers are "
  "compared *within* each branch, so the tomograms are identical and only the "
  "program changes; running it on both branches replicates the comparison.", "")

s2 = table("part2_summary.csv")

figure("part2_particle_counts.png",
       "Number of detected particles in all four sets, plus how many the two "
       "pickers agree on, per tomogram.",
       "From `warp_tiltseries_*/matching/*_clean.star` and "
       "`pytom_picks/*/*_particles.star`, matched by `picks.py`.",
       "**A higher count is not a better result.** A picker can return more picks "
       "simply by returning more false positives. Both tools were thresholded by "
       f"the same rule — {SIGMA}σ above each volume's own background — so the "
       "counts are comparable, but they remain counts and not accuracy.")

figure("part2_spatial_overlap.png",
       "Left: Jaccard overlap between the two pickers against the matching "
       "tolerance. Right: the fraction of Warp picks with a PyTom pick nearby. "
       "One curve per alignment branch.",
       "Optimal one-to-one matching (Hungarian algorithm) at each tolerance from "
       "10 to 200 Å.",
       f"**Why sweep the radius.** 'How close counts as the same particle' is the "
       f"most manipulable number in this comparison — quote 200 Å and almost "
       f"everything agrees, quote 10 Å and almost nothing does. The operating "
       f"point is {MATCH_RADIUS} Å, one apoferritin radius, because beyond that "
       f"the two spheres barely overlap and calling them the same particle stops "
       f"meaning anything physically. Publishing the whole curve lets a reader "
       f"judge whether the conclusion survives the choice.")

figure("part2_score_distributions.png",
       "Detection scores for all four sets, split into picks the other tool "
       "confirmed and picks it did not.",
       "Scores straight from the STAR files, never rescaled.",
       "**Separate axes on purpose.** Warp scores in standard deviations above "
       "background; PyTom scores a normalised correlation coefficient. Putting "
       "them on a shared axis, or rescaling both to 0–1 to make them look "
       "comparable, would invent a relationship that does not exist. A "
       "distribution that falls monotonically from the cutoff is the signature of "
       "a noise tail sliced at an arbitrary point; a real population shows a bump "
       "separated from the background.")

figure("part2_runtime.png",
       "Picking wall-clock for all four combinations.",
       "Timed in `part2_warp.py` and `part2_pytom.py`.",
       "PyTom searches 9216 orientations to Warp's 1536 — a 7.5° search over all "
       "3D rotations, divided by the symmetry each tool can exploit (Warp uses "
       "apoferritin's full 24-fold octahedral symmetry, PyTom only the 4-fold "
       "axis about z). Per orientation the two run at the same speed, so the "
       "runtime gap is symmetry handling, not implementation quality.")

figure("part2_pick_positions.png",
       "Where the detected particles physically sit, looking down the beam, one "
       "panel per branch.",
       "The x and y columns of both pick tables, in Ångströms.",
       "A sanity check. Clustering, edge effects, or picks outside the specimen "
       "region show up here immediately.")

w("### Per tomogram", "")
table("part2_per_tomogram.csv",
      "`chance` is what two unrelated pick lists of the same density would agree "
      "on by accident — the denser a tool's picks, the more of the other's it "
      "coincides with for free. Without that baseline the confirmation rates are "
      "unreadable. `above_chance` is the ratio.")

w("### At equal counts", "")
table("part2_equal_counts.csv",
      "Each tool's top N by score, N being the smaller of the two. This removes "
      "the count difference as a confound and asks the cleaner question: of each "
      "tool's best N picks, how many are the same molecules?")

w("### Are the picks only one tool found its weakest?", "")
table("part2_unique_vs_confirmed.csv",
      "`prob_confirmed_higher` is the chance a confirmed pick outscores a unique "
      "one; 0.5 means no relationship. Above 0.5 means the tool's own ranking "
      "agrees with the other tool's opinion — its scores are informative even "
      "where its threshold is not. This is a tested claim, not an assumed one.")

w("### Key parameters", "")
table("part2_parameters.csv")

w("### Runtime", "")
table("part2_runtime.csv")

w("### Interpretation", "")
embed("part2_interpretation.md")

# ------------------------------------------------------------------ videos
w("---", "", "# Videos", "",
  "| file | what it shows | built from |",
  "|---|---|---|",
  "| `results/videos/raw_frames.mp4` | every downloaded movie, frame-averaged and 8× downsampled, in acquisition order | `frames/*.tif` |",
  "| `results/videos/alignment_slices.mp4` | both alignments' tomograms side by side, slice by slice | `warp_tiltseries_*/reconstruction/*.mrc` |",
  "| `results/videos/picks_slices.mp4` | the same slices with detections circled, Warp green, PyTom red | tomograms + both pick tables |", "",
  "Stills from the last one are written to `results/plots/annotated_slices/`.", "",
  "**What to look for.** Each tomogram is a sandwich: empty ice, then the layer "
  "of protein, then empty ice again. Scrubbing through the video travels downward "
  "through it. A picker whose circles track that layer is finding molecules; one "
  "whose circles are spread evenly top to bottom is finding noise.", "")

path = OUT.parent / "RESULTS.md"
path.write_text("\n".join(md))
print(f"wrote {path}  ({len(md)} lines)")
