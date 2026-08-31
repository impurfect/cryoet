# Cryo-ET workflow comparison

Comparing two **tilt-series alignment** methods and two **particle-picking**
methods on the Warp tilt-series tutorial dataset (EMPIAR-10491, apoferritin),
with a reproducible analysis and a dashboard.

One script per job. No command-line arguments — every script runs as
`python <script>.py` and reads its settings from `config.py`.

---

# Part A — Setup

The three GPU programs (WarpTools, AreTomo2, PyTom) are **NVIDIA CUDA only**.
IMOD is Linux/macOS. There is no Windows build of WarpTools.

## Linux (Ubuntu 22.04 / 24.04, or Debian 12 / 13)

**1. NVIDIA driver.** On a GCP GPU instance:

```bash
sudo apt-get update
sudo apt-get install -y make gcc pciutils dkms cmake git "linux-headers-$(uname -r)"

curl -LO https://github.com/GoogleCloudPlatform/compute-gpu-installation/releases/download/cuda-installer-v1.9.1/cuda_installer.pyz
sudo mkdir -p /opt/google/cuda-installer && echo done | sudo tee /opt/google/cuda-installer/prerequisites
sudo python3 cuda_installer.pyz install_driver
nvidia-smi                                  # must list your GPU

sudo apt-mark hold "linux-image-$(uname -r)" "linux-headers-$(uname -r)"
```

The `prerequisites` line makes the installer skip its own kernel-upgrade step,
which picks kernel versions by sorting them as text and chooses one whose headers
do not exist. The `apt-mark hold` stops a later kernel upgrade from orphaning the
driver — a GPU VM reboots for host maintenance, and an unpinned kernel comes back
without a driver.

**2. System packages.**

```bash
sudo apt-get install -y build-essential wget curl unzip git python3-venv \
                        default-jre \
                        libgl1 libglu1-mesa libx11-6 libxext6 libxt6 libsm6 libice6
```

`default-jre` is not optional: IMOD's `etomo` is a Java program, and alignment
fails on every tilt series without it. The X and GL libraries are needed even for
IMOD's command-line tools.

**3. Conda environment** (one env holds everything):

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash && exec bash

conda create -n cryoet -y -c warpem -c pytorch -c conda-forge \
  warp=2.0.0 cupy "cuda-version=12.9" \
  numpy pandas scipy matplotlib mrcfile streamlit
conda activate cryoet
python -m pip install "pytom-match-pick[plotting]" -r requirements.txt
```

**4. IMOD and AreTomo2** (plain binaries, not conda packages):

```bash
# IMOD - the RHEL8 build is correct on Debian/Ubuntu; no .deb is published
cd /tmp && wget --no-check-certificate \
  https://bio3d.colorado.edu/imod/AMD64-RHEL5/imod_5.1.12_RHEL8-64_CUDA12.0.sh
sha256sum imod_5.1.12_RHEL8-64_CUDA12.0.sh
# expect 1cb30013c74f34a33313909cbaf293fb50fb07fa3cff71f2dec52d7b948c4da9
sudo sh imod_5.1.12_RHEL8-64_CUDA12.0.sh -yes
echo 'source /etc/profile.d/IMOD-linux.sh' >> ~/.bashrc && source /etc/profile.d/IMOD-linux.sh

# AreTomo2 - pick the CUDA 12 build to match the conda environment
mkdir -p ~/opt/aretomo2 && cd ~/opt/aretomo2
wget https://github.com/czimaginginstitute/AreTomo2/releases/download/v1.1.2/AreTomo2_1.1.2_03-27-2024.zip
unzip -o AreTomo2_1.1.2_03-27-2024.zip && chmod +x AreTomo2_*
ln -sf ~/opt/aretomo2/AreTomo2_1.1.2_Cuda121 ~/opt/aretomo2/AreTomo2
echo 'export PATH=$HOME/opt/aretomo2:$PATH' >> ~/.bashrc && source ~/.bashrc
```

AreTomo2 links against `libcudart.so.12` and `libcufft.so.11`, which live inside
the conda environment. It only starts with `cryoet` active — from a bare shell
you get `error while loading shared libraries` even though the install is fine.

The IMOD server sends an incomplete certificate chain, hence
`--no-check-certificate` plus a hash you can verify against.

**Ubuntu vs Debian:** identical. Same `apt`, same package names, same builds.
Only the GCP image family differs (`ubuntu-2204-lts` vs `debian-12`).

**RHEL / AlmaLinux / Rocky:** swap `apt-get install` for `dnf install`, use
`java-17-openjdk-headless`, `libGL libGLU libX11 libXext libXt mesa-libGL`, and
drop `python3-venv` (RHEL ships `venv` inside `python3`).

## Windows

**The processing cannot run natively on Windows.** WarpTools ships no Windows
binary, and AreTomo2 and PyTom are Linux/CUDA. Two options:

**Option 1 — WSL2 (Windows 11, or Windows 10 21H2+).** WSL2 passes an NVIDIA GPU
through to a Linux environment, and CUDA works inside it.

```powershell
wsl --install -d Ubuntu-22.04          # then reboot
```

Install the **Windows** NVIDIA driver on the host — not inside WSL, which would
break the passthrough. Then open the Ubuntu shell, confirm `nvidia-smi` works,
and follow the Linux instructions above from step 2. Skip step 1 entirely.

**Option 2 — analysis only.** `part1_analyze.py`, `part2_analyze.py`,
`part3_dashboard.py` and the video scripts need no GPU and no cryo-ET software.
Run the processing on a Linux machine, copy `results/` across, and work locally:

```powershell
py -3.11 -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python part1_analyze.py
streamlit run part3_dashboard.py
```

## Point the scripts at your data

```bash
export CRYOET_DATA=~/cryoet_data       # raw data and Warp's output folders
```

Everything else lives in `config.py`, which is the only file to edit.

---

# Part B — Prerequisites: what this is actually about

No background assumed. Read this before the assessment PDF and it will make
sense.

## The specimen — and what is *not* the target

**The ice is the background.** The sample is a solution containing millions of
copies of one protein, frozen so fast that the water becomes a glass rather than
crystals. The molecules are suspended in that glass, frozen in random positions
and orientations, like fruit set in jelly.

**The targets are the protein molecules** — apoferritin, an iron-storage protein.
In a tomogram slice they appear as small **rings or hexagons**, about 13 nm
across. That is what both picking programs are hunting for.

| word | what it means | what to look for in an image |
|---|---|---|
| **particle** | one apoferritin molecule | a small ring or hexagon |
| **ice** | the frozen water around them | the grainy grey background |
| **noise** | random speckle from the low electron dose | the fine "TV static" texture |
| **slice** / **z** | one flat layer through the 3D block | the `z=43` in a video caption |
| **contrast** | how much light-and-dark variation a slice has | many rings = high; flat grey = low |
| **pick** / **detection** | the program saying "a molecule is here" | one circle |
| **precision** | of the circles drawn, how many are on real molecules | are circles on rings, or on blank grey? |
| **recall** | of the real molecules, how many got circled | are there rings with no circle? |

## Why it is hard

Electrons destroy what they illuminate. So the dose is kept tiny — 2.64 electrons
per square Ångström per image — and **every raw image is almost pure noise**. You
genuinely cannot see a molecule in one. Everything downstream is about pulling
signal out from under the noise floor.

## From flat shadows to a 3D block

An electron micrograph is a **projection**: everything along the beam squashed
into one plane, like an X-ray of a hand. One projection cannot give you 3D; many
projections from different angles can. This is exactly a hospital CT scan, and
the mathematics is the same. The specimen is tilted step by step and photographed
at each angle.

- **Tilt series** — the set of 41 images taken at 41 angles, −40° to +40° in 2° steps
- **Tomogram** — the 3D volume reconstructed from them
- **Cryo-ET** — cryo-electron tomography, the whole technique

## `TS_1`, `TS_11`, `TS_17`, `TS_23`, `TS_32`

`TS` = Tilt Series. These are five **different patches of ice on the same grid**,
each containing different individual molecules of the same protein. Not the same
patch photographed five times.

The numbers are labels from the original experiment, not a count — the public
deposit holds 37 tilt series and the tutorial uses these five. They are five
independent samples, which is what makes "4 out of 5 tilt series agree" a
meaningful statement.

## 41 versus 78 — input versus output

These count entirely different things and are easy to confuse:

| number | what it is |
|---|---|
| **41** | camera angles — how many photographs were taken |
| **78** | slices — how many layers deep the reconstructed 3D block is |

A CT scanner takes X-rays from perhaps 40 angles, then the radiologist scrolls
through hundreds of cross-sections. The two numbers are unrelated.

Where 78 comes from — the ice thickness and how finely it is sampled, not the
number of images:

```
ice thickness   = 1000 raw pixels
raw pixel size  = 0.7894 Å
processing at   = 10 Å per voxel
1000 × 0.7894 / 10 = 79 slices
```

The full accounting:

```
INPUT    5 tilt series × 41 angles   = 205 raw images
OUTPUT   5 tomograms × 79 slices     = 395 slices per alignment branch
         × 2 branches                = 790 slices
```

Each volume is 347 × 474 nm across and 79 nm thick. So scrubbing through
`alignment_slices.mp4` means travelling **downward through the ice**: empty ice
at the top, the protein layer in the middle, empty ice at the bottom.

## What the whole exercise is for

You are not trying to photograph one molecule well — that is impossible at this
dose. The plan is to find **thousands of copies**, all frozen in different random
orientations, and average them. Random noise cancels; the shared shape
reinforces. That is how the tutorial reaches a 3 Å structure from data where no
individual molecule is visible.

Which is why particle picking matters, and why a picker that circles empty ice is
a real problem: feed those into the average and you are averaging static in with
signal.

## The four processing steps

1. **Motion correction** — the specimen creeps under the beam during each
   exposure; the frames are shifted back into register before averaging.
2. **Alignment** — the nominal stage angles are not accurate enough, so the true
   geometry is recovered from the images themselves. **Comparison 1.**
3. **Reconstruction** — back-project the aligned images into a 3D volume.
4. **Particle picking** — slide a known 3D shape (the *template*) through the
   volume and record where it fits. **Comparison 2.**

## The three questions

| | question | how |
|---|---|---|
| **1** | Does the alignment method change the science? | Process the same data twice, changing only the aligner, and carry both branches through to particle picking. |
| **2** | Does the particle picker change the science? | Run both pickers on the same tomograms with the same template and threshold. |
| **3** | Can someone else check the answer? | A dashboard driven entirely by the result tables. |

**Why both branches go all the way to picking.** IMOD and AreTomo each print an
"error" when they finish, and it is tempting to declare the smaller one the
winner. That comparison is invalid: IMOD reports the mean distance between
predicted and observed patch positions, AreTomo reports an error from a different
objective. Different quantities, different scales. So the alignments are judged on
outcomes measured after the branches rejoin — tomogram quality, and how
confidently molecules can be found.

---

# Part C — Running it

```bash
conda activate cryoet
python download.py          # ~3 GB, EMPIAR-10491 + EMD-15854; safe to re-run
python preprocess.py        # motion correction, CTF, tilt-series grouping
python part1_align.py       # both alignments + reconstruction
python part2_warp.py        # Warp template matching, both branches
python part2_pytom.py       # PyTom template matching, both branches
python part1_analyze.py     # Part 1 tables, plots, interpretation
python part2_analyze.py     # Part 2 tables, plots, interpretation
streamlit run part3_dashboard.py
```

Videos are optional and independent — run them whenever their inputs exist:

```bash
python raw_data_movie.py    # after download.py
python video_alignment.py   # after part1_align.py
python video_picks.py       # after both picking scripts
```

Rough timings on one NVIDIA L4: preprocessing 20 min, alignment 35 min, Warp
picking 3 min, PyTom picking 25 min, everything else seconds.

---

# Part D — What each script does

## `config.py`
Every path and parameter in one place: dataset identity, microscope values,
processing pixel size, the two branch folders, template, thresholds, GPU.
Nothing else hard-codes a number.

## `download.py`
Fetches the gain reference, five `.mdoc` metadata files, 205 movies and the
EMD-15854 template map.

The movie list comes from the `.mdoc` files, not a wildcard. The tutorial's own
script globs `*-11_*.tif`, which is ambiguous — the deposit holds both a `53-11`
and a `59-11` series — so that pattern quietly downloads eight tilt series where
you wanted five. Safe to re-run: `wget -N` skips what is already current and
never deletes.

## `preprocess.py`
The steps shared by both alignment branches, so they run once.

| command | argument | what it does |
|---|---|---|
| `create_settings` | `--folder_data frames` | where the raw movies are |
| | `--angpix 0.7894` | size of one camera pixel at the specimen, in Å |
| | `--gain_path gain_ref.mrc` | camera calibration image; every movie is divided by it to remove the sensor's fixed pattern |
| | `--gain_flip_y` | this camera's gain reference is stored upside down relative to the movies |
| | `--exposure 2.64` | electrons per Å² per image, used to track cumulative radiation damage |
| `fs_motion_and_ctf` | `--m_grid 1x1x3` | motion model resolution (X×Y×time). Drift is tracked over 3 time points but not across the image — there is too little signal per tilt to do more |
| | `--c_grid 2x2x1` | defocus model resolution. 2×2 across the image lets you check defocus varies sensibly along the tilt axis; constant in time |
| | `--c_range_max 7` | finest detail used for defocus fitting, in Å |
| | `--c_defocus_max 8` | largest defocus to consider, in µm |
| | `--c_use_sum` | fit defocus from the motion-corrected average rather than per-frame spectra — better when signal per frame is low |
| | `--out_averages` | write the drift-corrected average of each movie |
| `ts_import` | `--mdocs mdoc` | metadata telling Warp which movie belongs to which tilt series, at what angle |
| | `--tilt_exposure 2.64` | dose per tilt, for the cumulative-damage model |
| | `--min_intensity 0.3` | drop images darker than 30% of expected — usually the beam was blocked |
| | `--dont_invert` | sets the geometric handedness of the reconstruction. A property of this microscope; the tutorial establishes it for this dataset |
| `create_settings` | `--tomo_dimensions 4400x6000x1000` | size of the volume to reconstruct, in raw pixels. Z is the ice thickness |

## `part1_align.py` — Part 1

Aligns the same preprocessed data twice, then reconstructs both. Warp's
`--output_processing` sends each aligner's results to its own folder, so
everything before the branch is shared and everything after is identical.

**Method A — `ts_etomo_patches` (IMOD).** Cuts each image into squares and follows
them from tilt to tilt by cross-correlation, then fits one geometry explaining all
the tracks.

| argument | value | what it does |
|---|---|---|
| `--angpix` | 10 | work at 10 Å/pixel — fine enough to align, ~160× less data than full resolution |
| `--patch_size` | 500 | side of a tracking square, in Å. Large enough to contain trackable features, small enough that the specimen does not deform within one |
| `--initial_axis` | −85.6 | starting guess for the tilt-axis angle, refined from there |
| `--perdevice` | 2 | worker processes per GPU |
| `--output_processing` | `warp_tiltseries_etomo` | keeps this branch separate |

**Method B — `ts_aretomo` (AreTomo2).** Tracks nothing. Reconstructs a rough
volume, re-projects it to 2D, compares against the real images, corrects the
geometry, repeats.

| argument | value | what it does |
|---|---|---|
| `--alignz` | 800 | thickness (raw px) it assumes while iterating. Too small and it discards real specimen; too large and it fits noise |
| `--axis_iter` | 5 | tilt-axis refinement passes. **This dominates AreTomo's runtime** — five full passes before the final alignment. `1` is ~6× faster |
| `--min_fov` | 0 | keep every tilt. Above 0, tilts whose field of view overlaps too little get dropped |

**Then, per branch:**

| command | argument | what it does |
|---|---|---|
| `ts_defocus_hand` | `--check` | there is a sign ambiguity in how defocus varies across a tilted specimen; getting it backwards costs resolution. This measures it. Positive correlation = correct as-is |
| `ts_ctf` | `--range_high 7` | finest detail used for fitting, in Å |
| | `--defocus_max 8` | largest defocus considered, µm |
| `ts_reconstruct` | `--angpix 10` | build the 3D volume at 10 Å/voxel |

## `part2_warp.py` — Part 2, picker A

| command | argument | what it does |
|---|---|---|
| `ts_template_match` | `--tomo_angpix 10` | which tomograms to match against |
| | `--template_emdb 15854` | download the reference structure from the EMDB. Warp reads its voxel size from the map header automatically |
| | `--template_diameter 130` | size of the molecule, in Å. Also the minimum spacing between accepted peaks |
| | `--symmetry O` | octahedral. 24 rotations leave apoferritin unchanged, so only 1/24 of orientation space needs searching |
| | `--subdivisions 3` | angular step: 2 = 15°, **3 = 7.5°**, 4 = 3.75°. Finer is more accurate and much slower |
| | `--whiten` | spectral whitening — boost high frequencies so low-resolution power does not dominate the score |
| | `--perdevice 1` | one worker per GPU; template matching is memory-hungry |
| `threshold_picks` | `--minimum 3` | keep peaks at least 3 standard deviations above the volume's background |
| | `--in_suffix` / `--out_suffix` | which correlation volumes to read, and what to name the output |

Warp also applies `--max_missing_tilts 2` internally: positions not covered by at
least that many tilts are discarded. That is why it draws almost nothing in the
empty ice at the top and bottom of the volume.

## `part2_pytom.py` — Part 2, picker B

| command | argument | what it does |
|---|---|---|
| `pytom_create_template.py` | `-i emd_15854.map` | the same reference map Warp uses |
| | `--input-voxel-size-angstrom` | **read from the map header** (0.729 Å). Assuming a round 1.0 would rescale the template by 1.37× and it would match nothing |
| | `--output-voxel-size-angstrom 10` | resample to the tomogram's voxel size |
| | `--center` | centre the density by centre of mass |
| `pytom_create_mask.py` | `-b` | box size — **read from the template**, since PyTom requires the two to match exactly |
| | `-r 6.5` | mask radius in **voxels** (65 Å at 10 Å/voxel). Only what is inside contributes to the score, so surrounding noise does not dilute it |
| | `-s 1` | soft Gaussian edge, avoiding artefacts from a hard cutoff |
| `pytom_match_template.py` | `--particle-diameter 130` | molecule size, matching Warp |
| | `--angular-search 7.5` | set explicitly to equal Warp's 7.5°. PyTom's own default derives it from the diameter |
| | `--z-axis-rotational-symmetry 4` | PyTom supports symmetry only about z. Octahedral contains a 4-fold axis, and EMDB deposits O maps with it along z. So PyTom searches **9216** orientations to Warp's **1536** |
| | `--low-pass 30` | resolution limit, Å. 20 would be Nyquist, i.e. no filtering at all — and combined with whitening that amplifies the noisiest frequencies |
| | `--search-z z0 z1` | **restrict to the specimen slab**, found per tomogram from per-slice variance. Warp does this for itself; without it PyTom searches empty ice and reports thousands of detections there |
| | `--spectral-whitening` | matches Warp's `--whiten` |
| | `--random-phase-correction` | also match a phase-scrambled template, to calibrate the background. PyTom's own recommendation |
| | `--warp-xml-file` | hands PyTom Warp's tilt angles, dose and defocus, so both see identical metadata |
| `pytom_extract_candidates.py` | `--cut-off` | **computed as mean + 3σ of the correlation volume** — the same rule Warp applies. PyTom's automatic estimate instead fits a false-alarm model, which on this data chose a cutoff above every peak and extracted nothing |
| | `-n 3000` | ceiling on picks per tomogram, set high enough never to bind |
| | `--particle-diameter 130` | minimum spacing between accepted peaks |
| | `--tophat-filter` | *optional* (`PYTOM_TOPHAT` in config). Keeps only sharp peaks. This is tuning, not a correction — report it as "PyTom after tuning" alongside the default run |

## `picks.py`
Loads both pickers' STAR files into one table, in Ångströms, and provides the
matching routine.

**Why the conversion matters.** Warp writes coordinates *normalised to 0–1*
across the volume; PyTom writes *voxels*. Comparing the raw columns would place
every Warp particle within one voxel of the origin. Scores are left on their own
scales and never mixed.

**Matching** is solved as an assignment problem (Hungarian algorithm), not by
nearest neighbour. Nearest neighbour double-counts, and first-come-first-served
makes the answer depend on row order in the file.

## `part1_analyze.py` — Part 1 analysis
Writes `results/tables/part1_*.csv`, six plots and `part1_interpretation.md`.
The four metric groups the assessment names:

| output | metric | note |
|---|---|---|
| `part1_alignment_residuals.png` | **alignment residuals** | each method's own error. **Never compared between methods** — IMOD reports mean tracked-patch position error, AreTomo exposes no per-series equivalent, and they minimise different objectives |
| `part1_reconstruction_contrast.png` | **reconstruction quality** | spread of voxel values; a smeared tomogram tends to uniform grey |
| `part1_reconstruction_sharpness.png` | **reconstruction quality** | variance of the Laplacian, the standard blur measure |
| `part1_particles_found.png` | downstream effect | particles per series, plus yield vs cutoff |
| `part1_peak_scores.png` | downstream effect | how strongly molecules stood out — the decisive metric |
| `part1_runtime.png` | **runtime** | wall clock for 5 tilt series |

Both reconstruction measures locate the specimen by per-slice variance rather
than assuming it sits mid-volume: the two aligners place it about 80 Å apart in
z, so a fixed slice would compare one branch's sample against the other's ice.

Comparisons are **paired** across the five tilt series. With five pairs a Wilcoxon
test cannot reach p<0.05 even when every pair agrees — 0.0625 is its floor — so
the evidence is consistency of direction, reported as "4/5 series".

## `part2_analyze.py` — Part 2 analysis
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

Three deliberate choices:

- **Sweeps the matching tolerance** from 10 to 200 Å. "How close counts as the
  same particle" is the most manipulable number in the comparison — quote 200 Å
  and everything agrees, quote 10 Å and nothing does. The operating point is
  65 Å, one apoferritin radius, with the whole curve published beside it.
- **Reports a chance baseline.** The denser a tool's picks, the more of the
  other's it coincides with for free. Without `chance` and `above_chance` the
  confirmation rates are unreadable.
- **Compares at equal counts too** — each tool's top N by score, N the smaller of
  the two — so the count difference stops driving the overlap statistic.

## `part3_dashboard.py` — Part 3
Streamlit page: alignment metrics, picking metrics, conclusions, videos. Every
number is read from `results/`, so the page always reflects the current run.

## `plotstyle.py`
Shared matplotlib settings and colours, so the plotting scripts hold only their
own logic.

## `video.py` and the three video scripts
`video.py` holds shared frame helpers. Each script writes one MP4 at 10 fps (set
`FPS` at the top):

- **`raw_data_movie.py`** — every downloaded movie, frame-averaged and 8×
  downsampled, in acquisition order. What the microscope actually recorded.
- **`video_alignment.py`** — both alignments' tomograms side by side, slice by
  slice. Same data, same reconstruction: any difference in sharpness is the
  aligner.
- **`video_picks.py`** — the same slices with each picker's detections circled,
  Warp green, PyTom red. A pick is drawn when its centre is within one particle
  radius of the slice. Also writes every 10th slice as a PNG into
  `results/plots/annotated_slices/` (`SAVE_EVERY` at the top).

---

# Part E — Outputs

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
├── part1_interpretation.md       written from the tables
├── part2_interpretation.md
└── runtime_*.json
```

---

# Part F — Notes on the data

- **EMPIAR-10491**, five tilt series of **apoferritin** — a hollow protein shell
  130 Å across with octahedral symmetry, the standard cryo-EM test specimen.
- **EMD-15854** is that same apoferritin, used as the template by both pickers.
  It is not a ribosome, and the dataset is not EMPIAR-10164 (immature HIV
  particles) — both are common mix-ups, and diameter, mask size, symmetry and
  match tolerance all follow from getting them right.
- Peaks reach roughly 4–5σ, close to the ~4.3σ ceiling noise alone would produce
  in a volume this size. A meaningful fraction of picks will be false positives.
  That is what five tilt series of a small particle at 10 Å/px gives you, and it
  is why the comparison rests on **paired differences between methods on
  identical data**, where a shared false-positive rate largely cancels.
- The definitive test — subtomogram averaging the agreed picks and comparing the
  resolution reached — is not done here. It needs RELION and M, and is the
  obvious next step.
