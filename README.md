# Cryo-ET workflow comparison

Comparing two **tilt-series alignment** methods and two **particle-picking**
methods on the Warp tilt-series tutorial dataset (EMPIAR-10491, apoferritin),
with a reproducible analysis and a dashboard.

One script per job. No command-line arguments — every script runs as
`python <script>.py` and reads its settings from `config.py`.

---

## The three tasks

| | question | how |
|---|---|---|
| **1** | Does the alignment method change the science? | Process the same data twice, changing only the aligner, and carry both branches through to particle picking. |
| **2** | Does the particle picker change the science? | Run both pickers on the same tomograms with the same template and threshold. |
| **3** | Can someone else check the answer? | A dashboard driven entirely by the result tables. |

**Why both branches go all the way to picking.** IMOD and AreTomo each print an
"error" when they finish, and it is tempting to declare the smaller one the
winner. That comparison is invalid: IMOD reports the mean distance between
predicted and observed patch positions, AreTomo reports an error from a different
projection-matching objective. They are different quantities on different scales.
So the alignments are judged on outcomes measured after the branches rejoin —
tomogram quality, and how confidently molecules can be found.

---

## Setup

Needs a Linux machine with an NVIDIA GPU. WarpTools, AreTomo2, IMOD and PyTom
all live in one conda environment; see the install notes in
`../cryoet_comparison/README.md` §2 if you are starting from scratch.

```bash
conda activate cryoet
pip install -r requirements.txt

export CRYOET_DATA=~/cryoet_data      # where raw data and Warp's output live
```

Everything in `config.py` is a plain constant — dataset, pixel size, dose,
template, thresholds, GPU. It is the only file to edit.

`analyze` and `dashboard` scripts need no GPU and no cryo-ET software: copy
`results/` to any machine and run them there.

---

## Running it, in order

```bash
python download.py          # ~3 GB, EMPIAR-10491 + EMD-15854; safe to re-run
python preprocess.py        # motion correction, CTF, tilt-series grouping
python part1_align.py       # both alignments + reconstruction
python part2_warp.py        # Warp template matching, both branches
python part2_pytom.py       # PyTom template matching, both branches
python part1_analyze.py     # Task 1 tables, plots, interpretation
python part2_analyze.py     # Task 2 tables, plots, interpretation
streamlit run part3_dashboard.py
```

Videos are optional and independent — run them whenever the inputs exist:

```bash
python raw_data_movie.py    # after download.py
python video_alignment.py   # after part1_align.py
python video_picks.py       # after both picking scripts
```

---

## What each script does

### `config.py`
Every path and parameter, in one place. Dataset identity, microscope values
(0.7894 Å/px, 2.64 e⁻/Å² per tilt), processing pixel size (10 Å/px), the two
alignment branch folders, template (EMD-15854, 130 Å, octahedral), the 3σ peak
cutoff and the 65 Å match radius. Nothing else hard-codes a number.

### `download.py`
Fetches the gain reference, the five `.mdoc` metadata files, the movies and the
EMD-15854 template map.

Safe to re-run: `wget -N` skips anything already on disk with an up-to-date
timestamp, so pointing it at an existing `cryoet_data` re-verifies rather than
re-downloads. It never deletes.

The movie list comes from the `.mdoc` files rather than a wildcard. The tutorial's
own script globs `*-11_*.tif`, which is ambiguous — this deposit contains both a
`53-11` and a `59-11` series — so that pattern quietly downloads eight tilt series
where you wanted five. The `.mdoc` files name exactly the 205 movies needed.

### `preprocess.py`
The steps shared by both alignment branches, so they run once:

- **Motion correction** — the sample creeps under the beam during each exposure;
  the frames are shifted back into register before averaging.
- **CTF estimation** — the microscope applies an oscillating transfer function
  that depends on defocus; it has to be measured from the power spectrum before
  it can be corrected.
- **Tilt-series grouping** — reads the `.mdoc` files and works out which movies
  belong to which tilt series, at what angle and accumulated dose.

### `part1_align.py` — Task 1
Aligns the same preprocessed data twice, then reconstructs both. Warp's
`--output_processing` sends each aligner's results to its own folder, so
everything before the branch point is shared and everything after is identical.

- **IMOD patch tracking** cuts each image into 500 Å squares and follows them
  between tilts by cross-correlation.
- **AreTomo2** never tracks anything: it reconstructs, re-projects to 2D,
  compares with the real images, corrects, and iterates.

Then per branch: defocus handedness check, tilt-series CTF estimation, and
reconstruction at 10 Å/px. Runtimes are written to
`results/runtime_alignment.json`.

### `part2_warp.py` — Task 2
Warp's built-in 3D template matching against EMD-15854, on **both** branches.
Correlates the template at every position and orientation on a 7.5° grid,
exploiting apoferritin's octahedral symmetry, then extracts peaks above 3σ.

### `part2_pytom.py` — Task 2
PyTom template matching on the same tomograms, from the same EMD-15854 map.

Two details that matter:

- The template's input voxel size is read from the map header. EMD-15854 is
  0.729 Å/voxel — assuming a round 1.0 rescales the template by 1.37× and it then
  matches nothing at all.
- The extraction cutoff is set to **3σ above each correlation volume's
  background**, the same rule Warp applies. PyTom's default instead fits a
  false-alarm model, which on this data chose a cutoff above every peak in the
  volume. Two pickers judged by two different rules cannot be compared on how
  many particles they find.

### `picks.py`
Loads both pickers' STAR files into one table, in Ångströms, and provides the
matching routine.

**Why the conversion matters.** Warp writes coordinates *normalised to 0–1*
across the volume; PyTom writes *voxels*. Comparing the raw columns would place
every Warp particle within one voxel of the origin. Scores are left on their own
scales and never mixed: Warp's is standard deviations above background, PyTom's
is a correlation coefficient.

**Matching** is solved as an assignment problem (Hungarian algorithm), not by
nearest neighbour. Nearest neighbour double-counts, and first-come-first-served
makes the answer depend on the order rows appear in the file.

### `part1_analyze.py` — Part 1 analysis
Writes `results/tables/part1_*.csv`, six plots, and `part1_interpretation.md`.
The four metric groups the assessment names:

| output | metric | note |
|---|---|---|
| `part1_alignment_residuals.png` | **alignment residuals** | each method's own error. **Never compared between methods** — IMOD reports mean tracked-patch position error, AreTomo exposes no per-series equivalent, and the two minimise different objectives |
| `part1_reconstruction_contrast.png` | **reconstruction quality** | spread of voxel values; a smeared tomogram tends to uniform grey |
| `part1_reconstruction_sharpness.png` | **reconstruction quality** | variance of the Laplacian, the standard blur measure |
| `part1_particles_found.png` | downstream effect | particles found per series, plus yield vs cutoff |
| `part1_peak_scores.png` | downstream effect | how strongly molecules stood out — the decisive metric |
| `part1_runtime.png` | **runtime** | wall clock for 5 tilt series |

Both reconstruction measures locate the specimen by per-slice variance rather
than assuming it sits mid-volume: the two aligners place it about 80 Å apart in
z, so a fixed slice would compare one branch's sample against the other's ice.

Comparisons are **paired** across the five tilt series. With five pairs a
Wilcoxon test cannot reach p<0.05 even when every pair agrees — 0.0625 is its
floor — so the evidence is consistency of direction, reported as "4/5 series".

### `part2_analyze.py` — Part 2 analysis
Two pickers × two alignment branches = **four pick sets**, and every metric
reports all four. The pickers are compared *within* each branch, which is the
controlled experiment — identical tomograms, only the program changes — and
running it on both branches replicates that comparison.

| output | metric |
|---|---|
| `part2_particle_counts.png` | **number of detected particles**, all four sets, plus how many agree |
| `part2_spatial_overlap.png` | **spatial overlap**, per branch, against matching tolerance |
| `part2_score_distributions.png` | **detection score distributions**, 2 branches × 2 pickers |
| `part2_runtime.png` | **runtime**, all four |
| `part2_pick_positions.png` | where the picks physically sit, per branch |
| `part2_parameters.csv` | **key parameters**, side by side |

Three things it does deliberately:

- **Sweeps the matching tolerance** from 10 to 200 Å. "How close counts as the
  same particle" is the most manipulable number in the comparison — quote 200 Å
  and everything agrees, quote 10 Å and nothing does. The operating point is
  65 Å, one apoferritin radius, and the whole curve is published beside it.
- **Reports a chance baseline.** The denser a tool's picks, the more of the
  other's it coincides with for free. Without `chance` and `above_chance` the
  confirmation rates are unreadable.
- **Compares at equal counts too** — each tool's top N by score, N the smaller
  of the two — so the count difference stops driving the overlap statistic.

### `part3_dashboard.py` — Task 3
Streamlit page: alignment metrics, picking metrics, conclusions, and the videos.
Every number is read from `results/` — nothing is written into the file, so the
page always reflects the current run.

### `plotstyle.py`
Shared matplotlib settings and the colour map, so the plotting scripts contain
only their own logic.

### `video.py`, `raw_data_movie.py`, `video_alignment.py`, `video_picks.py`
`video.py` holds the shared frame helpers. The three scripts each write one MP4
at 10 fps (set `FPS` at the top of each):

- **`raw_data_movie.py`** — every downloaded movie, frame-averaged and
  8× downsampled, in acquisition order. What the microscope actually recorded.
- **`video_alignment.py`** — the two alignments' tomograms side by side, slice by
  slice. Same data, same reconstruction: any difference in sharpness is the
  aligner.
- **`video_picks.py`** — the same slices with each picker's detections circled,
  Warp in green, PyTom in red. A pick is drawn on a slice when its centre is
  within one particle radius, so each particle appears across the slices it
  spans. As well as the MP4 it writes every 10th slice as a PNG into
  `results/plots/annotated_slices/` (change `SAVE_EVERY` at the top), which is
  what you want for a report or a slide.

---

## Outputs

```
results/
├── tables/                       every measurement, as CSV
│   ├── part1_summary.csv         paired comparison of the two alignments
│   ├── part1_per_series.csv      the raw numbers behind it
│   ├── part1_alignment_residuals.csv
│   ├── part2_summary.csv         headline picking numbers, per branch
│   ├── part2_per_tomogram.csv    counts and overlap, all four sets
│   ├── part2_radius_sweep.csv    agreement vs matching tolerance, per branch
│   ├── part2_equal_counts.csv    the same comparison with counts equalised
│   ├── part2_runtime.csv         all four runtimes
│   ├── part2_parameters.csv      what each picker was told to do
│   └── part2_unique_vs_confirmed.csv
├── plots/                        11 PNGs
│   └── annotated_slices/         picks circled on tomogram slices
├── videos/                       3 MP4s
├── task1_interpretation.md       written from the tables
├── task2_interpretation.md
└── runtime_*.json
```

---

## Notes on the data

- **EMPIAR-10491**, five tilt series of **apoferritin** — a hollow protein shell
  130 Å across with octahedral symmetry, the standard cryo-EM test specimen.
- **EMD-15854** is that same apoferritin, used as the template by both pickers.
  It is not a ribosome, and the dataset is not EMPIAR-10164 (immature HIV
  particles) — both are common mix-ups, and diameter, mask size, symmetry and
  match tolerance all follow from getting them right.
- Peaks in this data reach roughly 4–5σ, close to the ~4.3σ ceiling that noise
  alone would produce in a volume this size. A meaningful fraction of the picks
  will be false positives. That is what five tilt series of a small particle at
  10 Å/px gives you, and it is why the comparison rests on **paired differences
  between methods on identical data**, where a shared false-positive rate largely
  cancels, rather than on any absolute count.
- The definitive test — subtomogram averaging the agreed picks and comparing the
  resolution reached — is not done here. It needs RELION and M, and is the
  obvious next step.
