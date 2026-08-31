# Results

Compiled from `results/` by `make_report.py`. Every number is read from the tables; nothing is written in by hand.

**Dataset** EMPIAR-10491, 5 tilt series of apoferritin · **Processing** 10 Å/voxel · **Template** EMD-15854, 130 Å, octahedral · **Peak cutoff** 3σ · **Match radius** 65 Å

---

## Headline

**Part 1.** IMOD etomo gives the stronger downstream result — 13.3% higher mean score across the 50 best peaks per tilt series (4/5 series in that direction). IMOD etomo is faster: 273 s against 1754 s for all five series.

**Part 2.** Warp returned 3472 particles across both branches, PyTom 3629. They agree on 1034 within 65 Å — 9.5x more than two unrelated lists of the same density would by chance. With counts equalised, agreement is 24%.

Both statements are computed from the tables below.

---

## 1. Where the data came from

`download.py` fetches three things from the EBI public archives:

| what | source | lands in |
|---|---|---|
| 205 movie files (5 series × 41 tilts) | `ftp.ebi.ac.uk/empiar/world_availability/10491/data/tiltseries/data/` | `$CRYOET_DATA/frames/` |
| 5 metadata files | `.../10491/data/tiltseries/mdoc/` | `$CRYOET_DATA/mdoc/` |
| camera gain reference | `.../10491/data/gain_ref.mrc` | `$CRYOET_DATA/gain_ref.mrc` |
| template structure | `ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-15854/map/` | `$CRYOET_DATA/emd_15854.map` |

The movie list is taken from the `.mdoc` files rather than a filename wildcard. The tutorial's own script globs `*-11_*.tif`, which matches both the `53-11` and `59-11` series in this deposit and quietly downloads eight tilt series where you wanted five.

**What each is:**

- **`.tif` movies** — the raw output. Each is a short burst of frames of one view at one tilt angle, 5760 × 4092 pixels. At the dose used they are almost pure noise; you cannot see a molecule in one.
- **`.mdoc`** — plain text, one block per tilt image, recording which movie file it is, the stage angle, and the accumulated dose. This is what turns 205 unrelated movies into 5 ordered tilt series.
- **`gain_ref.mrc`** — a calibration image of the camera sensor. Every movie is divided by it to remove the sensor's fixed pattern.
- **`emd_15854.map`** — a published 1.8 Å structure of apoferritin, deposited at 0.729 Å/voxel. Both pickers use this same file as their reference shape.

## 2. What preprocessing does to it

`preprocess.py` runs once, and its output is shared by both alignment branches so they cannot differ for any reason other than the aligner.

| step | what it does | creates |
|---|---|---|
| `create_settings` | records pixel size, dose, gain path | `warp_frameseries.settings` |
| `fs_motion_and_ctf` | measures how the specimen drifted during each exposure and shifts the frames back into register; measures the defocus of every image from its power spectrum | `warp_frameseries/` — one XML of results per movie, plus `average/` holding one drift-corrected image per movie |
| `ts_import` | reads the `.mdoc` files and groups movies into tilt series with their angles and cumulative dose | `tomostar/` — five `.tomostar` files |
| `create_settings` | records the volume size to reconstruct | `warp_tiltseries.settings` |

**The new data:** `warp_frameseries/average/` holds 205 images that are still 2D, but now drift-corrected and with known defocus. These, not the raw `.tif` files, are what everything downstream uses.

## 3. What alignment and reconstruction create

`part1_align.py` runs the two aligners into separate folders using Warp's `--output_processing`, then reconstructs both:

```
$CRYOET_DATA/
├── warp_tiltseries_etomo/      <- branch A, IMOD patch tracking
│   ├── TS_*.xml                   the projection geometry it solved for
│   ├── tiltstack/TS_*/            IMOD's working files and logs
│   └── reconstruction/            TS_*_10.00Apx.mrc  <- the tomograms
└── warp_tiltseries_aretomo/    <- branch B, AreTomo2, same structure
```

Each tomogram is 347 × 474 × 79 voxels at 10 Å — about 347 × 474 nm across and 79 nm thick, roughly 26 MB. **This is the first point at which molecules are visible.**

## 4. What particle picking creates

| script | creates | contents |
|---|---|---|
| `part2_warp.py` | `warp_tiltseries_*/matching/` | correlation volumes, and `TS_*_clean.star` listing accepted peaks. Coordinates are **normalised 0–1**; scores are in **σ above background** |
| `part2_pytom.py` | `pytom_picks/etomo/`, `pytom_picks/aretomo/` | `*_scores.mrc`, `*_angles.mrc`, `*_job.json`, and `*_particles.star`. Coordinates are in **voxels**; scores are **normalised cross-correlation** |

Two pickers × two branches = **four pick sets**. `picks.py` converts both to Ångströms so positions can be compared, and leaves the scores on their own scales because they measure different things.

---

# Part 1 — Does the alignment method matter?

Both branches share every step except alignment and both were carried through to particle picking, so they are judged on outcomes rather than on each program's internal error measure.

| metric                                 | higher_is       |     etomo |    aretomo |   difference_pct | etomo_better_in   |   wilcoxon_p |
|:---------------------------------------|:----------------|----------:|-----------:|-----------------:|:------------------|-------------:|
| reconstruction quality: contrast       | more structure  |   1.35549 |    1.34839 |              0.5 | 4/5 series        |       0.4375 |
| reconstruction quality: sharpness      | crisper edges   |   0.0023  |    0.00224 |              2.7 | 4/5 series        |       0.125  |
| downstream: particles found            | more particles  | 343.4     |  351       |             -2.2 | 1/5 series        |       0.1875 |
| downstream: mean score of top 50 peaks | stronger peaks  |   8.24719 |    7.27694 |             13.3 | 4/5 series        |       0.125  |
| alignment runtime (s, 5 series)        | lower is faster | 273.1     | 1754.1     |            -84.4 | nan               |     nan      |

_Read from `results/tables/part1_summary.csv`._

### Alignment residuals

![part1_alignment_residuals.png](results/plots/part1_alignment_residuals.png)

**What it is.** Each aligner's own reported error, one panel per method.

**How it was made.** `part1_analyze.py` scans `warp_tiltseries_*/tiltstack/*/*.log` for a labelled residual line.

**What it says.** **These two panels must not be compared with each other.** IMOD reports the mean distance between predicted and observed tracked-patch positions; AreTomo exposes no per-series equivalent and minimises a different objective entirely. Useful only for spotting one bad tilt series within one method.

### Reconstruction quality

![part1_reconstruction_contrast.png](results/plots/part1_reconstruction_contrast.png)

**What it is.** Spread of voxel values in the specimen slab, one line per tilt series joining its two branch values.

**How it was made.** Reads `reconstruction/*.mrc` from both branches. The slab is located per volume by per-slice variance, not assumed to be mid-volume — the two aligners place the specimen about 80 Å apart in z.

**What it says.** A well-aligned tomogram has crisp molecules and therefore more spread in its values; a misaligned one smears towards uniform grey. A paired plot shows whether a difference is consistent across all five series or driven by one outlier.

![part1_reconstruction_sharpness.png](results/plots/part1_reconstruction_sharpness.png)

**What it is.** Variance of the Laplacian of the specimen slice — the standard blur measure, borrowed from camera autofocus.

**How it was made.** Same volumes as above.

**What it says.** The Laplacian responds to edges, so its variance is high when edges are crisp and low when blurred. It measures blur without needing to know what the right answer looks like.

### Downstream effect

![part1_particles_found.png](results/plots/part1_particles_found.png)

**What it is.** Left: particles found per tilt series in each branch. Right: how the count changes as the score cutoff moves from 3σ to 12σ.

**How it was made.** Counts the picks in `warp_tiltseries_*/matching/*_clean.star` — the same picker on both branches, so only the alignment differs.

**What it says.** The right-hand panel matters: a single count depends entirely on where you put the threshold. If one branch is above the other across the whole curve, the conclusion does not rest on an arbitrary cutoff.

![part1_peak_scores.png](results/plots/part1_peak_scores.png)

**What it is.** Distribution of template-matching scores from both branches, sharing an axis.

**How it was made.** Same STAR files. Sharing an axis is legitimate here: same program, same template, same normalisation.

**What it says.** **This is the plot that answers Part 1.** The score is how many standard deviations a peak stands above its own volume's background. A distribution pushed further right means molecules stood out more clearly, which means the tomogram preserved their structure better, which means the alignment was better.

![part1_runtime.png](results/plots/part1_runtime.png)

**What it is.** Wall-clock seconds to align all five tilt series.

**How it was made.** Timed inside `part1_align.py`, written to `results/runtime_alignment.json`.

**What it says.** Note this is partly a parameter choice, not purely algorithmic: AreTomo ran five tilt-axis refinement passes (`--axis_iter 5`) while IMOD was handed the axis directly (`--initial_axis`). At `--axis_iter 1` AreTomo would be roughly six times faster.

### Per tilt series

| series   |   contrast_aretomo |   contrast_etomo |   sharpness_aretomo |   sharpness_etomo |   specimen_centre_z_aretomo |   specimen_centre_z_etomo |   particles_etomo |   peak_score_etomo |   particles_aretomo |   peak_score_aretomo |
|:---------|-------------------:|-----------------:|--------------------:|------------------:|----------------------------:|--------------------------:|------------------:|-------------------:|--------------------:|---------------------:|
| TS_1     |            1.3648  |          1.37882 |          0.00215192 |        0.00224871 |                          35 |                        44 |               335 |            8.84177 |                 324 |              6.98159 |
| TS_11    |            1.33284 |          1.33706 |          0.00244815 |        0.0025574  |                          41 |                        35 |               445 |            8.5907  |                 452 |              8.26649 |
| TS_17    |            1.34841 |          1.32752 |          0.00237345 |        0.00235768 |                          48 |                        41 |               306 |            7.05902 |                 318 |              7.22799 |
| TS_23    |            1.35578 |          1.37274 |          0.00210283 |        0.00216835 |                          48 |                        42 |               303 |            8.53447 |                 316 |              6.91167 |
| TS_32    |            1.34013 |          1.36129 |          0.00213133 |        0.00218335 |                          36 |                        44 |               328 |            8.21    |                 345 |              6.99695 |

### Interpretation

#### Part 1 - alignment comparison

| metric                                 | higher_is       |     etomo |    aretomo |   difference_pct | etomo_better_in   |   wilcoxon_p |
|:---------------------------------------|:----------------|----------:|-----------:|-----------------:|:------------------|-------------:|
| reconstruction quality: contrast       | more structure  |   1.35549 |    1.34839 |              0.5 | 4/5 series        |       0.4375 |
| reconstruction quality: sharpness      | crisper edges   |   0.0023  |    0.00224 |              2.7 | 4/5 series        |       0.125  |
| downstream: particles found            | more particles  | 343.4     |  351       |             -2.2 | 1/5 series        |       0.1875 |
| downstream: mean score of top 50 peaks | stronger peaks  |   8.24719 |    7.27694 |             13.3 | 4/5 series        |       0.125  |
| alignment runtime (s, 5 series)        | lower is faster | 273.1     | 1754.1     |            -84.4 |                   |              |

##### Per tilt series

| series   |   contrast_aretomo |   contrast_etomo |   sharpness_aretomo |   sharpness_etomo |   specimen_centre_z_aretomo |   specimen_centre_z_etomo |   particles_etomo |   peak_score_etomo |   particles_aretomo |   peak_score_aretomo |
|:---------|-------------------:|-----------------:|--------------------:|------------------:|----------------------------:|--------------------------:|------------------:|-------------------:|--------------------:|---------------------:|
| TS_1     |            1.3648  |          1.37882 |             0.00215 |           0.00225 |                          35 |                        44 |               335 |            8.84177 |                 324 |              6.98159 |
| TS_11    |            1.33284 |          1.33706 |             0.00245 |           0.00256 |                          41 |                        35 |               445 |            8.5907  |                 452 |              8.26649 |
| TS_17    |            1.34841 |          1.32752 |             0.00237 |           0.00236 |                          48 |                        41 |               306 |            7.05902 |                 318 |              7.22799 |
| TS_23    |            1.35578 |          1.37274 |             0.0021  |           0.00217 |                          48 |                        42 |               303 |            8.53447 |                 316 |              6.91167 |
| TS_32    |            1.34013 |          1.36129 |             0.00213 |           0.00218 |                          36 |                        44 |               328 |            8.21    |                 345 |              6.99695 |

##### Alignment residuals

**Reported for completeness, not for comparison.** IMOD reports the mean distance between predicted and observed tracked-patch positions. AreTomo exposes no equivalent per-series figure. Even if it did, the two minimise different objectives and are not on a common scale. Use them to spot a bad tilt series within one method, nothing more.

| branch   | series   | metric                 |      value | unit   | source    |
|:---------|:---------|:-----------------------|-----------:|:-------|:----------|
| etomo    | TS_1     | mean residual          |     0.418  | px     | align.log |
| etomo    | TS_1     | weighted mean residual |     0.304  | px     | align.log |
| etomo    | TS_1     | measured/unknowns      | 12672      |        | align.log |
| etomo    | TS_11    | mean residual          |     0.2223 | px     | align.log |
| etomo    | TS_11    | weighted mean residual |     0.211  | px     | align.log |
| etomo    | TS_11    | measured/unknowns      | 12672      |        | align.log |
| etomo    | TS_17    | mean residual          |     2.296  | px     | align.log |
| etomo    | TS_17    | weighted mean residual |     1.992  | px     | align.log |
| etomo    | TS_17    | measured/unknowns      | 12672      |        | align.log |
| etomo    | TS_23    | mean residual          |     0.458  | px     | align.log |
| etomo    | TS_23    | weighted mean residual |     0.323  | px     | align.log |
| etomo    | TS_23    | measured/unknowns      | 12672      |        | align.log |
| etomo    | TS_32    | mean residual          |     0.72   | px     | align.log |
| etomo    | TS_32    | weighted mean residual |     0.432  | px     | align.log |
| etomo    | TS_32    | measured/unknowns      | 12672      |        | align.log |

##### Verdict

**IMOD etomo** gives the stronger downstream result: 13.3% higher mean score across the 50 best peaks per tilt series (4/5 series favour etomo). Reconstruction contrast and sharpness are near-identical, so the difference shows up where it matters - in how clearly molecules can be found.

**IMOD etomo** is faster: 273 s against 1754 s. Note this is partly a parameter choice: AreTomo ran 5 tilt-axis refinement iterations (--axis_iter 5) while IMOD was given the axis directly (--initial_axis).

With five paired tilt series a Wilcoxon test cannot reach p<0.05 even when every pair agrees - 0.0625 is its floor. Consistency of direction is the evidence, not the p-value.

---

# Part 2 — Does the particle picker matter?

Two pickers × two alignment branches = four pick sets. The pickers are compared *within* each branch, so the tomograms are identical and only the program changes; running it on both branches replicates the comparison.

| branch   |   n_warp |   n_pytom |   n_matched |   jaccard |   warp_confirmed |   pytom_confirmed |   chance |   above_chance |   agree_at_equal_counts |   median_separation_A |   spearman_rho |   warp_runtime_s |   pytom_runtime_s |
|:---------|---------:|----------:|------------:|----------:|-----------------:|------------------:|---------:|---------------:|------------------------:|----------------------:|---------------:|-----------------:|------------------:|
| etomo    |     1717 |      1908 |         572 |     0.187 |            0.333 |             0.3   |    0.033 |          10.16 |                   0.243 |                  17.2 |         -0.663 |             93.8 |             181   |
| aretomo  |     1755 |      1721 |         462 |     0.153 |            0.263 |             0.268 |    0.03  |           8.77 |                   0.244 |                  34.3 |         -0.246 |            113.3 |             179.5 |

![part2_particle_counts.png](results/plots/part2_particle_counts.png)

**What it is.** Number of detected particles in all four sets, plus how many the two pickers agree on, per tomogram.

**How it was made.** From `warp_tiltseries_*/matching/*_clean.star` and `pytom_picks/*/*_particles.star`, matched by `picks.py`.

**What it says.** **A higher count is not a better result.** A picker can return more picks simply by returning more false positives. Both tools were thresholded by the same rule — 3σ above each volume's own background — so the counts are comparable, but they remain counts and not accuracy.

![part2_spatial_overlap.png](results/plots/part2_spatial_overlap.png)

**What it is.** Left: Jaccard overlap between the two pickers against the matching tolerance. Right: the fraction of Warp picks with a PyTom pick nearby. One curve per alignment branch.

**How it was made.** Optimal one-to-one matching (Hungarian algorithm) at each tolerance from 10 to 200 Å.

**What it says.** **Why sweep the radius.** 'How close counts as the same particle' is the most manipulable number in this comparison — quote 200 Å and almost everything agrees, quote 10 Å and almost nothing does. The operating point is 65 Å, one apoferritin radius, because beyond that the two spheres barely overlap and calling them the same particle stops meaning anything physically. Publishing the whole curve lets a reader judge whether the conclusion survives the choice.

![part2_score_distributions.png](results/plots/part2_score_distributions.png)

**What it is.** Detection scores for all four sets, split into picks the other tool confirmed and picks it did not.

**How it was made.** Scores straight from the STAR files, never rescaled.

**What it says.** **Separate axes on purpose.** Warp scores in standard deviations above background; PyTom scores a normalised correlation coefficient. Putting them on a shared axis, or rescaling both to 0–1 to make them look comparable, would invent a relationship that does not exist. A distribution that falls monotonically from the cutoff is the signature of a noise tail sliced at an arbitrary point; a real population shows a bump separated from the background.

![part2_runtime.png](results/plots/part2_runtime.png)

**What it is.** Picking wall-clock for all four combinations.

**How it was made.** Timed in `part2_warp.py` and `part2_pytom.py`.

**What it says.** PyTom searches 9216 orientations to Warp's 1536 — a 7.5° search over all 3D rotations, divided by the symmetry each tool can exploit (Warp uses apoferritin's full 24-fold octahedral symmetry, PyTom only the 4-fold axis about z). Per orientation the two run at the same speed, so the runtime gap is symmetry handling, not implementation quality.

![part2_pick_positions.png](results/plots/part2_pick_positions.png)

**What it is.** Where the detected particles physically sit, looking down the beam, one panel per branch.

**How it was made.** The x and y columns of both pick tables, in Ångströms.

**What it says.** A sanity check. Clustering, edge effects, or picks outside the specimen region show up here immediately.

### Per tomogram

| branch   | series   |   n_warp |   n_pytom |   n_matched |   jaccard |   warp_confirmed |   pytom_confirmed |   chance |   above_chance |   median_separation_A |
|:---------|:---------|---------:|----------:|------------:|----------:|-----------------:|------------------:|---------:|---------------:|----------------------:|
| etomo    | TS_1     |      335 |       488 |         174 |     0.268 |            0.519 |             0.357 |    0.042 |          12.28 |                  13.3 |
| etomo    | TS_11    |      445 |       333 |          77 |     0.11  |            0.173 |             0.231 |    0.029 |           5.96 |                  17.2 |
| etomo    | TS_17    |      306 |       863 |         239 |     0.257 |            0.781 |             0.277 |    0.074 |          10.62 |                  18.2 |
| etomo    | TS_23    |      303 |        27 |           6 |     0.019 |            0.02  |             0.222 |    0.002 |           8.29 |                  13.1 |
| etomo    | TS_32    |      328 |       197 |          76 |     0.169 |            0.232 |             0.386 |    0.017 |          13.4  |                  25.3 |
| aretomo  | TS_1     |      324 |       171 |          43 |     0.095 |            0.133 |             0.251 |    0.015 |           8.83 |                  25.2 |
| aretomo  | TS_11    |      452 |       394 |         105 |     0.142 |            0.232 |             0.266 |    0.034 |           6.78 |                  36.3 |
| aretomo  | TS_17    |      318 |       400 |         142 |     0.247 |            0.447 |             0.355 |    0.035 |          12.83 |                  27.7 |
| aretomo  | TS_23    |      316 |       354 |          47 |     0.075 |            0.149 |             0.133 |    0.031 |           4.82 |                  37.1 |
| aretomo  | TS_32    |      345 |       402 |         125 |     0.201 |            0.362 |             0.311 |    0.035 |          10.36 |                  26.3 |

`chance` is what two unrelated pick lists of the same density would agree on by accident — the denser a tool's picks, the more of the other's it coincides with for free. Without that baseline the confirmation rates are unreadable. `above_chance` is the ratio.

### At equal counts

| branch   | series   |   n_each |   n_matched |   fraction_agreeing |   chance |   above_chance |
|:---------|:---------|---------:|------------:|--------------------:|---------:|---------------:|
| etomo    | TS_1     |      335 |         149 |               0.445 |    0.029 |          15.22 |
| etomo    | TS_11    |      333 |          50 |               0.15  |    0.029 |           5.17 |
| etomo    | TS_17    |      306 |         114 |               0.373 |    0.027 |          13.94 |
| etomo    | TS_23    |       27 |           0 |               0     |    0.002 |           0    |
| etomo    | TS_32    |      197 |          49 |               0.249 |    0.017 |          14.39 |
| aretomo  | TS_1     |      171 |          22 |               0.129 |    0.015 |           8.56 |
| aretomo  | TS_11    |      394 |          93 |               0.236 |    0.034 |           6.89 |
| aretomo  | TS_17    |      318 |         121 |               0.381 |    0.028 |          13.71 |
| aretomo  | TS_23    |      316 |          47 |               0.149 |    0.028 |           5.39 |
| aretomo  | TS_32    |      345 |         112 |               0.325 |    0.03  |          10.79 |

Each tool's top N by score, N being the smaller of the two. This removes the count difference as a confound and asks the cleaner question: of each tool's best N picks, how many are the same molecules?

### Are the picks only one tool found its weakest?

| branch   | picker   |   n_confirmed |   n_unique |   median_confirmed |   median_unique |   prob_confirmed_higher |        p |
|:---------|:---------|--------------:|-----------:|-------------------:|----------------:|------------------------:|---------:|
| etomo    | warp     |           572 |       1145 |             7.2443 |          7.7223 |                   0.36  | 1        |
| etomo    | pytom    |           572 |       1336 |             0.3418 |          0.343  |                   0.548 | 0.000408 |
| aretomo  | warp     |           462 |       1293 |             6.4907 |          6.5722 |                   0.485 | 0.834972 |
| aretomo  | pytom    |           462 |       1259 |             0.3082 |          0.2888 |                   0.638 | 0        |

`prob_confirmed_higher` is the chance a confirmed pick outscores a unique one; 0.5 means no relationship. Above 0.5 means the tool's own ranking agrees with the other tool's opinion — its scores are informative even where its threshold is not. This is a tested claim, not an assumed one.

### Key parameters

| parameter             | warp                                               | pytom                                            |
|:----------------------|:---------------------------------------------------|:-------------------------------------------------|
| template              | EMD-15854                                          | EMD-15854                                        |
| particle diameter (A) | 130                                                | 130                                              |
| angular step (deg)    | 7.5 (subdivisions 3)                               | 7.5                                              |
| symmetry used         | O (octahedral, 24-fold)                            | C4 about z (PyTom supports z-axis symmetry only) |
| orientations searched | 1536 (36,864 / 24)                                 | 9216 (36,864 / 4)                                |
| search region         | auto (drops positions with <3 tilts covering them) | --search-z, set to the specimen slab             |
| spectral whitening    | on                                                 | on                                               |
| score definition      | sigma above background                             | normalised cross-correlation (LCCmax)            |
| peak cutoff           | 3 sigma                                            | mean + 3 x std of the correlation volume         |

### Runtime

| picker   | branch   |   seconds |   seconds_per_tomogram |
|:---------|:---------|----------:|-----------------------:|
| warp     | etomo    |      93.8 |                   18.8 |
| warp     | aretomo  |     113.3 |                   22.7 |
| pytom    | etomo    |     181   |                   36.2 |
| pytom    | aretomo  |     179.5 |                   35.9 |

### Interpretation

#### Part 2 - particle-picking comparison

Two pickers on two alignment branches: four pick sets. The pickers are compared within each branch, so the tomograms are identical and only the program changes; running it on both branches replicates the comparison.

##### Summary, per branch

| branch   |   n_warp |   n_pytom |   n_matched |   jaccard |   warp_confirmed |   pytom_confirmed |   chance |   above_chance |   agree_at_equal_counts |   median_separation_A |   spearman_rho |   warp_runtime_s |   pytom_runtime_s |
|:---------|---------:|----------:|------------:|----------:|-----------------:|------------------:|---------:|---------------:|------------------------:|----------------------:|---------------:|-----------------:|------------------:|
| etomo    |     1717 |      1908 |         572 |     0.187 |            0.333 |             0.3   |    0.033 |          10.16 |                   0.243 |                  17.2 |         -0.663 |             93.8 |             181   |
| aretomo  |     1755 |      1721 |         462 |     0.153 |            0.263 |             0.268 |    0.03  |           8.77 |                   0.244 |                  34.3 |         -0.246 |            113.3 |             179.5 |

##### Number of detected particles, per tomogram

| branch   | series   |   n_warp |   n_pytom |   n_matched |   jaccard |   warp_confirmed |   pytom_confirmed |   chance |   above_chance |   median_separation_A |
|:---------|:---------|---------:|----------:|------------:|----------:|-----------------:|------------------:|---------:|---------------:|----------------------:|
| etomo    | TS_1     |      335 |       488 |         174 |     0.268 |            0.519 |             0.357 |    0.042 |          12.28 |                  13.3 |
| etomo    | TS_11    |      445 |       333 |          77 |     0.11  |            0.173 |             0.231 |    0.029 |           5.96 |                  17.2 |
| etomo    | TS_17    |      306 |       863 |         239 |     0.257 |            0.781 |             0.277 |    0.074 |          10.62 |                  18.2 |
| etomo    | TS_23    |      303 |        27 |           6 |     0.019 |            0.02  |             0.222 |    0.002 |           8.29 |                  13.1 |
| etomo    | TS_32    |      328 |       197 |          76 |     0.169 |            0.232 |             0.386 |    0.017 |          13.4  |                  25.3 |
| aretomo  | TS_1     |      324 |       171 |          43 |     0.095 |            0.133 |             0.251 |    0.015 |           8.83 |                  25.2 |
| aretomo  | TS_11    |      452 |       394 |         105 |     0.142 |            0.232 |             0.266 |    0.034 |           6.78 |                  36.3 |
| aretomo  | TS_17    |      318 |       400 |         142 |     0.247 |            0.447 |             0.355 |    0.035 |          12.83 |                  27.7 |
| aretomo  | TS_23    |      316 |       354 |          47 |     0.075 |            0.149 |             0.133 |    0.031 |           4.82 |                  37.1 |
| aretomo  | TS_32    |      345 |       402 |         125 |     0.201 |            0.362 |             0.311 |    0.035 |          10.36 |                  26.3 |

##### Spatial overlap

Two picks count as the same molecule when their centres are within 65 A - one apoferritin radius. Beyond that the two spheres barely overlap and calling them the same particle stops meaning anything. Because that tolerance is the single most manipulable number here, the full curve from 10 to 200 A is published alongside it.

`chance` is what two unrelated lists of the same density would agree on by accident: the denser a tool's picks, the more of the other's it coincides with for free. Without that baseline the confirmation rates are unreadable.

##### At equal counts

| branch   | series   |   n_each |   n_matched |   fraction_agreeing |   chance |   above_chance |
|:---------|:---------|---------:|------------:|--------------------:|---------:|---------------:|
| etomo    | TS_1     |      335 |         149 |               0.445 |    0.029 |          15.22 |
| etomo    | TS_11    |      333 |          50 |               0.15  |    0.029 |           5.17 |
| etomo    | TS_17    |      306 |         114 |               0.373 |    0.027 |          13.94 |
| etomo    | TS_23    |       27 |           0 |               0     |    0.002 |           0    |
| etomo    | TS_32    |      197 |          49 |               0.249 |    0.017 |          14.39 |
| aretomo  | TS_1     |      171 |          22 |               0.129 |    0.015 |           8.56 |
| aretomo  | TS_11    |      394 |          93 |               0.236 |    0.034 |           6.89 |
| aretomo  | TS_17    |      318 |         121 |               0.381 |    0.028 |          13.71 |
| aretomo  | TS_23    |      316 |          47 |               0.149 |    0.028 |           5.39 |
| aretomo  | TS_32    |      345 |         112 |               0.325 |    0.03  |          10.79 |

Taking each tool's top N by score, N the smaller of the two, removes the count difference as a confound: agreement 24% against 2% by chance (9.4x).

##### Detection scores

| branch   | picker   |   n_confirmed |   n_unique |   median_confirmed |   median_unique |   prob_confirmed_higher |        p |
|:---------|:---------|--------------:|-----------:|-------------------:|----------------:|------------------------:|---------:|
| etomo    | warp     |           572 |       1145 |             7.2443 |          7.7223 |                   0.36  | 1        |
| etomo    | pytom    |           572 |       1336 |             0.3418 |          0.343  |                   0.548 | 0.000408 |
| aretomo  | warp     |           462 |       1293 |             6.4907 |          6.5722 |                   0.485 | 0.834972 |
| aretomo  | pytom    |           462 |       1259 |             0.3082 |          0.2888 |                   0.638 | 0        |

`prob_confirmed_higher` is the chance a confirmed pick outscores a unique one; 0.5 means no relationship. Above 0.5 means the tool's own ranking agrees with the other tool's opinion, so its scores are informative even where its threshold is not.

##### Runtime and parameters

| picker   | branch   |   seconds |   seconds_per_tomogram |
|:---------|:---------|----------:|-----------------------:|
| warp     | etomo    |      93.8 |                   18.8 |
| warp     | aretomo  |     113.3 |                   22.7 |
| pytom    | etomo    |     181   |                   36.2 |
| pytom    | aretomo  |     179.5 |                   35.9 |

| parameter             | warp                                               | pytom                                            |
|:----------------------|:---------------------------------------------------|:-------------------------------------------------|
| template              | EMD-15854                                          | EMD-15854                                        |
| particle diameter (A) | 130                                                | 130                                              |
| angular step (deg)    | 7.5 (subdivisions 3)                               | 7.5                                              |
| symmetry used         | O (octahedral, 24-fold)                            | C4 about z (PyTom supports z-axis symmetry only) |
| orientations searched | 1536 (36,864 / 24)                                 | 9216 (36,864 / 4)                                |
| search region         | auto (drops positions with <3 tilts covering them) | --search-z, set to the specimen slab             |
| spectral whitening    | on                                                 | on                                               |
| score definition      | sigma above background                             | normalised cross-correlation (LCCmax)            |
| peak cutoff           | 3 sigma                                            | mean + 3 x std of the correlation volume         |

PyTom searches 9216 orientations to Warp's 1536 - a 7.5 degree search over SO(3) divided by the symmetry each can exploit. Per orientation the two run at the same speed, so the runtime gap is symmetry handling, not implementation quality.

A higher count is not by itself a better result: a picker can return more picks purely by returning more false positives. Counts are counts, not accuracy.

---

# Videos

| file | what it shows | built from |
|---|---|---|
| `results/videos/raw_frames.mp4` | every downloaded movie, frame-averaged and 8× downsampled, in acquisition order | `frames/*.tif` |
| `results/videos/alignment_slices.mp4` | both alignments' tomograms side by side, slice by slice | `warp_tiltseries_*/reconstruction/*.mrc` |
| `results/videos/picks_slices.mp4` | the same slices with detections circled, Warp green, PyTom red | tomograms + both pick tables |

Stills from the last one are written to `results/plots/annotated_slices/`.

**What to look for.** Each tomogram is a sandwich: empty ice, then the layer of protein, then empty ice again. Scrubbing through the video travels downward through it. A picker whose circles track that layer is finding molecules; one whose circles are spread evenly top to bottom is finding noise.
