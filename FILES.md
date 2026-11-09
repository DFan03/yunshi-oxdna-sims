# Data drive map (FILES)

## 1. What is on this drive

This drive holds simulation data (RNA/ and DNA/) and wet-lab AFM images (WetLab/). RNA/ has 14 simulation runs in 8 design/size groups, plus starting structures, umbrella-sampling and forward-flux (ffs) folders. DNA/ has 7 runs in 6 groups, umbrella-sampling folders and nanotweezer figures. WetLab/AFM has 98 files from 9 measurement days. This page describes the folder layout only, not the scripts (scripts are not described here: see repo docs/rna.md and repo docs/dna.md). Generated 2026-11-09 12:46 from the manifests.

| Top folder | Files | Total size | Source of the number |
|---|---|---|---|
| RNA | 1962 | 5.7 GB | RNA/MANIFEST.csv |
| DNA (placed) | 3826 | 12.3 GB | DNA/MANIFEST.csv, without unsorted _raw |
| DNA/_raw (still arriving) | 2 | 783.0 MB | DNA/MANIFEST.csv, role raw_unsorted (transfer running; do not use) |
| WetLab | 98 | 347.0 MB | file sizes on disk (no manifest) |

## 2. Quick find: where is X?

| Question | Path | Section |
|---|---|---|
| Production trajectory of RNA design D, run N | `RNA/simulations/<D>/<monomer\|8mer>/run_NN/03_production_MD/`; example: `RNA/simulations/2H_standard/8mer/run_01/03_production_MD/trajectory.dat.zst` | 4 |
| Production trajectory of a DNA design | `DNA/simulations/<D>/<monomer\|8mer\|3mer\|design>/run_NN/03_production_MD/` | 4 |
| Umbrella window trajectories and set-up | `DNA/umbrella/<source>/` and `RNA/umbrella/<source>/` (window_NN/equil, prod) | 4 |
| Umbrella free-energy (WHAM) profiles | `DNA/umbrella/umbrella_DNA/umbrella`, `DNA/umbrella/umbrella_DNA/umbrella_500confs`, `DNA/umbrella/umbrella_DNA/umbrella_k1/umbrella`, `DNA/umbrella/umbrella_DNA_opening/umbrella`, `DNA/umbrella/umbrella_K1IVT8mer/umbrella`, `DNA/umbrella/umbrella_K1IVT8mer/umbrella_utils` (+6 more) | 4 |
| Starting structures (RNA) | `RNA/starting_structures/` | 4 |
| Starting structures (nanotweezer) | `DNA/simulations/nanotw_starting_structures/` | 4 |
| Forward-flux sampling (ffs) | `RNA/ffs/` | 8 |
| Nanotweezer figures and angle csv | `DNA/nanotw_figures/` | 8 |
| Wet-lab AFM of date D | `WetLab/AFM/YYYYMMDD/`; days present: 2026-08-26, 2026-08-28, 2026-09-01, 2026-09-02, 2026-09-03, 2026-09-07, 2026-09-08, 2026-09-09, 2026-09-10 | 5 |
| What a file type is | table in section 7 | 7 |
| What is missing or uncertain | section 6 | 6 |

## 3. Layout of one run

```
<RNA|DNA>/simulations/<design>/<monomer|8mer|3mer|design>/run_NN/
    01_relaxation_MC/   - Monte Carlo relaxation: first, gentle clean-up of the start structure (inputs, last structure).
    02_relaxation_MD/   - Short molecular-dynamics relaxation after MC; its last structure feeds production.
    03_production_MD/   - The main run: trajectory.dat.zst, last_conf.dat, energy.dat, input.txt, structure.top, log.
    analysis/   - Analysis outputs (csv, plots, json, mean structures, aligned trajectory) made from the trajectory.
    scripts_used/   - Copies of the helper scripts that were in the run folder (not explained here).
    starting_structure/   - Root .dat/.top the run was built from (DNA nanotweezer runs).
    other/   - Files that did not fit anywhere else (nested inside a stage folder).
    README.md, ORIGINAL_NAMES.tsv
```

Folders found in the manifests (number of runs that have them): 01_relaxation_MC 21, 02_relaxation_MD 21, 03_production_MD 21, analysis 21, other 6, scripts_used 21, starting_structure 3. Run numbering: run_01 is the newest main production run; run_02 is older. Total runs listed: 21.

## 4. Simulation runs and umbrella sampling

### RNA designs (monomer / 8mer)

| Design | Size | Run | Date range | Status | Files | Size on disk | Trajectory file | Supersedes / notes |
|---|---|---|---|---|---|---|---|---|
| 2H_standard | 8mer | run_01 | 2026-03-29 .. 2026-08-12 | production trajectory present | 57 | 983.4 MB | `trajectory.dat.zst` (541.5 MB) | supersedes: none; superseded by: none |
| 2H_standard | monomer | run_01 | 2026-03-18 .. 2026-03-30 | production trajectory present | 47 | 187.9 MB | `trajectory.dat.zst` (67.3 MB) | supersedes: none; superseded by: none |
| K1_IVT_curved | 8mer | run_01 | 2026-03-29 .. 2026-08-12 | production trajectory present | 62 | 893.5 MB | `trajectory.dat.zst` (497.4 MB) | supersedes: run_02; superseded by: none |
| K1_IVT_curved | 8mer | run_02 | 2026-03-10 .. 2026-03-12 | no production trajectory (input only) | 14 | 7.6 MB | none | supersedes: none; superseded by: run_01 |
| K1_IVT_curved | monomer | run_01 | 2026-03-18 .. 2026-05-27 | production trajectory present | 53 | 173.4 MB | `trajectory.dat.zst` (62.6 MB) | supersedes: run_02; superseded by: none |
| K1_IVT_curved | monomer | run_02 | 2026-03-09 .. 2026-03-10 | production trajectory present | 35 | 65.0 MB | `trajectory.dat.zst` (61.0 MB) | supersedes: none; superseded by: run_01 |
| K1_curved | 8mer | run_01 | 2026-03-29 .. 2026-03-31 | production trajectory present | 46 | 903.0 MB | `trajectory.dat.zst` (494.1 MB) | supersedes: run_02; superseded by: none; design name unconfirmed |
| K1_curved | 8mer | run_02 | 2026-03-10 .. 2026-03-11 | no production trajectory (input only) | 16 | 10.2 MB | none | supersedes: none; superseded by: run_01; design name unconfirmed |
| K1_curved | monomer | run_01 | 2026-03-18 .. 2026-05-21 | production trajectory present | 51 | 174.1 MB | `trajectory.dat.zst` (62.4 MB) | supersedes: run_02; superseded by: none; design name unconfirmed |
| K1_curved | monomer | run_02 | 2026-03-09 .. 2026-03-10 | production trajectory present | 38 | 142.4 MB | `trajectory.dat.zst` (77.5 MB) | supersedes: none; superseded by: run_01; design name unconfirmed |
| T3_IVT_wavy | monomer | run_01 | 2026-03-18 .. 2026-04-22 | production trajectory present | 49 | 188.9 MB | `trajectory.dat.zst` (66.2 MB) | supersedes: run_02; superseded by: none |
| T3_IVT_wavy | monomer | run_02 | 2026-03-09 .. 2026-03-11 | production trajectory present | 38 | 67.3 MB | `trajectory.dat.zst` (64.7 MB) | supersedes: none; superseded by: run_01 |
| T3_wavy | monomer | run_01 | 2026-03-18 .. 2026-04-22 | production trajectory present | 51 | 188.1 MB | `trajectory.dat.zst` (65.8 MB) | supersedes: run_02; superseded by: none; design name unconfirmed |
| T3_wavy | monomer | run_02 | 2026-03-09 .. 2026-03-11 | production trajectory present | 36 | 67.4 MB | `trajectory.dat.zst` (64.8 MB) | supersedes: none; superseded by: run_01; design name unconfirmed |

### DNA designs

| Design | Size | Run | Date range | Status | Files | Size on disk | Trajectory file | Supersedes / notes |
|---|---|---|---|---|---|---|---|---|
| T3_IVT_wavy | 3mer | run_01 | 2026-05-04 .. 2026-05-04 | relax only (no production MD in the source) | 26 | 63.1 MB | none | run_01 (only one source per DNA design here) |
| T3_IVT_wavy | 8mer | run_01 | 2026-05-05 .. 2026-11-09 | production trajectory present | 70 | 5.5 GB | `trajectory.dat.zst` present | run_01 (only one source per DNA design here) |
| T3_IVT_wavy | monomer | run_01 | 2026-04-13 .. 2026-11-09 | production trajectory present (03_production_MD/trajectory.dat.zst) | 51 | 693.2 MB | `trajectory.dat.zst` (674.7 MB) | run_01 (only one source per DNA design here) |
| chol_fluoro_nanotw | design | run_01 | 2026-06-17 .. 2026-06-17 | relax only (no production MD in the source) | 23 | 129.0 MB | none | run_01 (only one source per DNA design here) |
| chol_fluoro_nanotw | design | run_02 | 2026-06-17 .. 2026-06-17 | relax only (no production MD in the source) | 20 | 67.2 MB | none | run_02 (newest = run_01 rule: the db copy in run_01 has the later mtimes, this wsl copy is older -> run_02 |
| nanotw_closed_patch4 | design | run_01 | 2026-04-13 .. 2026-09-14 | production MD run present; trajectory in 03_production_MD | 33 | 439.8 MB | `trajectory.dat.zst` (337.8 MB) | run_01 |
| nanotw_open_patch4 | design | run_01 | 2026-03-19 .. 2026-09-14 | production MD run present; trajectory in 03_production_MD | 34 | 3.1 GB | `trajectory.dat.zst` (1.5 GB); `trajectory_complete.dat.zst` (1.5 GB) | run_01 |

Date range = earliest to latest file date in the run (copy times of .zst files excluded where the README says so). Size 'design' = multiplicity not specified.

| open_patch4 trajectory | Size | Last timestep | Frames |
|---|---|---|---|
| `trajectory.dat.zst` | 1.5 GB | 452,000,000 | 904 |
| `trajectory_complete.dat.zst` | 1.5 GB | 451,500,000 | 903 |

### Umbrella sampling

| Source | Where | What it is | Files | Size | Windows | Reaction coordinate range | Window trajectories kept | Large files left out? |
|---|---|---|---|---|---|---|---|---|
| nanotw_umbrella_sampling | DNA/umbrella/nanotw_umbrella_sampling | umbrella-sampling setup/test for the nanotweezer | 15 | 7.9 MB | 5 | 40.0-64.0 | equil 0, prod 0 (of 5 windows) | no file >50 MB was present |
| nanotw_ustest | DNA/umbrella/nanotw_ustest | umbrella-sampling setup/test for the nanotweezer | 415 | 279.4 MB | 25 | 27-61 | equil 33, prod 1 (of 25 windows) | no file >50 MB was present |
| umbrella_DNA | DNA/umbrella/umbrella_DNA | umbrella-sampling window runs, OUTSIDE the numbered run_NN layout | 1339 | 698.1 MB | not recorded | not recorded | equil 86, prod 10 | yes: files >50 MB left out |
| umbrella_DNA_opening | DNA/umbrella/umbrella_DNA_opening | umbrella-sampling window runs, OUTSIDE the numbered run_NN layout | 435 | 284.9 MB | 25 | 60-9 | equil 25, prod 0 (of 25 windows) | yes: files >50 MB left out |
| umbrella_K1IVT8mer | DNA/umbrella/umbrella_K1IVT8mer | umbrella-sampling window runs, OUTSIDE the numbered run_NN layout | 411 | 45.0 MB | 25 | 7-12 | equil 25, prod 0 (of 25 windows) | yes: files >50 MB left out |
| umbrella_T3IVT | DNA/umbrella/umbrella_T3IVT | umbrella-sampling window runs, OUTSIDE the numbered run_NN layout | 467 | 897.6 MB | 25 | 7.2-12 | equil 25, prod 25 (of 25 windows) | yes: files >50 MB left out |
| umbrella_T3IVT8mer | DNA/umbrella/umbrella_T3IVT8mer | umbrella-sampling window runs, OUTSIDE the numbered run_NN layout | 433 | 47.1 MB | 25 | 7.3-12 | equil 25, prod 0 (of 25 windows) | yes: files >50 MB left out |
| umbrella_K1IVT8mer | RNA/umbrella/umbrella_K1IVT8mer | umbrella-sampling run (no README_SORT.md) | 406 | 44.8 MB | 25 | 7-12 | equil 25, prod 0 (of 25 windows) | unknown |
| us_k1 | RNA/umbrella/us_k1 | umbrella-sampling run (no README_SORT.md) | 332 | 683.3 MB | 20 | 7.5-10.2 | equil 20, prod 20 (of 20 windows) | unknown |

Window trajectories larger than 50 MB were NOT copied (transfer rule in each README_SORT.md); the exact names of the skipped files are not recorded. Setup files, inputs, energies, small trajectories and profiles are included.

## 5. Wet lab: AFM images

| Date | Files | Size | File types | Measurement series | Empty folders |
|---|---|---|---|---|---|
| 2026-08-26 | 1 | 850.0 kB | spm 1 | 1 | New folder |
| 2026-08-28 | 10 | 14.9 MB | spm 10 | 1 | - |
| 2026-09-01 | 3 | 7.9 MB | spm 3 | 1 | - |
| 2026-09-02 | 15 | 28.5 MB | spm 15 | 2 | New folder |
| 2026-09-03 | 27 | 95.8 MB | spm 27 | 4 | - |
| 2026-09-07 | 12 | 26.0 MB | spm 9, svg 3 | 3 | - |
| 2026-09-08 | 9 | 65.4 MB | spm 9 | 1 | - |
| 2026-09-09 | 3 | 7.3 MB | spm 3 | 1 | - |
| 2026-09-10 | 18 | 100.2 MB | spm 18 | 3 | - |

Sample table: the columns are copied from the file name only. Time tokens (like 40min, 2H) are copied as written; their exact meaning is not interpreted here. 'unparsed' means the file name does not say.

| Date | Series (file name stem) | Strand | Time tokens | Stage words | Files | Types |
|---|---|---|---|---|---|---|
| 2026-08-26 | `S1C+C_ReOpen_halfwashed` | S1 | unparsed | ReOpen_halfwashed | 1 | spm 1 |
| 2026-08-28 | `S1C+C(no_toggle)` | S1 | unparsed | no_toggle | 10 | spm 10 |
| 2026-09-01 | `S2C+C_Air` | S2 | unparsed | unparsed | 3 | spm 3 |
| 2026-09-02 | `S1C+C_S2C+C_assmble_together_2h` | S1, S2 | 2h | assmble_together_2h | 14 | spm 14 |
| 2026-09-02 | `S2C+C_Air` | S2 | unparsed | unparsed | 1 | spm 1 |
| 2026-09-03 | `S1C+C_1h_S2C+C_1h_48c2HX2_RTovernight` | S1, S2 | 1h, 1h, 2H | unparsed | 7 | spm 7 |
| 2026-09-03 | `S1C+C_1h_S2C+C_1h_48c2HX3_RTovernight` | S1, S2 | 1h, 1h, 2H | unparsed | 3 | spm 3 |
| 2026-09-03 | `S1C+C_1h_S2C+C_1h_48c2H_RTovernight` | S1, S2 | 1h, 1h, 2H | unparsed | 9 | spm 9 |
| 2026-09-03 | `S1C+C_ReOpen2H` | S1 | 2H | ReOpen2H | 8 | spm 8 |
| 2026-09-07 | `S1C+C40min` | S1 | 40min | unparsed | 3 | spm 3 |
| 2026-09-07 | `S1C+C40min_S2C+C40min_3X2Hassmble_reopen1H` | S1, S2 | 40min, 40min, 2H, 1H | assmble_reopen1H | 1 | spm 1 |
| 2026-09-07 | `S1C+C_S2C+C_assemble2X_reopen1H` | S1, S2 | 1H | assemble2X_reopen1H | 8 | spm 5, svg 3 |
| 2026-09-08 | `S1C+C40min_S2C+C40min_3X2Hassemble_2HreOpen` | S1, S2 | 40min, 40min, 2H, 2H | assemble_2HreOpen | 9 | spm 9 |
| 2026-09-09 | `S1+C` | S1 | unparsed | unparsed | 3 | spm 3 |
| 2026-09-10 | `S1+C_S2+C_3Hassmble` | S1, S2 | 3H | assmble | 5 | spm 5 |
| 2026-09-10 | `S1+C_S2+C_3hassmble_reclose1h` | S1, S2 | 3h, 1h | assmble_reclose1h | 8 | spm 8 |
| 2026-09-10 | `S1+C_reclose3h` | S1 | 3h | reclose3h | 5 | spm 5 |

No sidecar text files (notes) exist next to the images, so the table relies on file names only.

## 6. Read this before using

> **Read this before using**
>
> - **.zst files (37 here) must be decompressed first:** `zstd -d file.zst` (keeps the .zst). Details: [READ_ME_zst.txt](RNA/READ_ME_zst.txt).
> - Decompressed trajectories are several times larger; the disk is nearly full. Decompress only what you need.
> - A file is only 'there' if it is listed in the MANIFEST.csv of RNA/ or DNA/.

### Deliberately left out or not placed

- Umbrella window trajectories larger than 50 MB (sources: umbrella_DNA, umbrella_DNA_opening, umbrella_K1IVT8mer, umbrella_T3IVT, umbrella_T3IVT8mer). Exact file names not recorded.
- RNA `tests/` (rna_ys/tests): not on this drive yet; will arrive later (RNA/SORT_NOTES.md).
- Unfinished or stale transfer copies kept in DNA/_raw (do not use): `_raw/T3IVT/T3IVT_md/T3IVT_trj.dat` (390.6 MB), `_raw/T3IVT_8mer/T3IVT_8mer_md/T3IVT_8mer_md.dat` (392.4 MB)
- Duplicates: `*.backup` files and the close_path4 `analysis/trj.dat` (reported identical by md5 to existing trajectories) were not placed (DNA/SORT_NOTES.md).
- `rna_utils` and `dna_utils` helper folders: not part of the data map; scripts live in the git repo.
- `align.dat(.zst)` aligned trajectories sit in `analysis/aligned/` (8 found); they are not in the original spec, so check before reuse.

### Uncertain or unparsed items

- **WetLab:** 12 of 51 sample-table cells (strand/time/stage) are 'unparsed': the file names do not say.
- **WetLab:** Meaning of 'C+C', 'X2/X3', '48c' in AFM file names is not recorded; tokens copied literally.
- **Model:** interaction_type=RNA2 (oxRNA model) appears in DNA inputs of: T3_IVT_wavy. Not checked whether intended.
- **T3_IVT_wavy 8mer:** Meaning of 'Mlink' (T3IVT_8mer_Mlink.dat/.top start files) is not recorded.
- **nanotw_open_patch4:** trajectory.dat.zst ends at t=452,000,000 (904 frames); trajectory_complete.dat.zst ends at t=451,500,000 (903 frames). They differ; the name 'complete' is not evidence; neither reaches the planned steps.
- **Umbrella:** us_k1 / umbrella_k1 are classed as RNA only because of the name (K1 RNA design).
- **chol_fluoro_nanotw:** run_02 (wsl copy) is older than run_01 (db copy); 9 files differ in size; relax only, no production.
- **T3_IVT_wavy 8mer:** trajectory.dat.zst: present
- **RNA designs:** Design names unconfirmed: K1_curved, T3_wavy.
- **RNA runs:** No production trajectory (only inputs): K1_IVT_curved/8mer/run_02, K1_curved/8mer/run_02.
- **Starting structures:** Underlying PDB/model of the start structure is 'unknown' in 19 run READMEs.
- **Umbrella:** Names of skipped (>50 MB) window trajectories are not recorded; per-window gaps are inferred.
- **Umbrella:** umbrella_K1IVT8mer exists in both RNA/ and DNA/ (DNA copy has 5 extra files, RNA copy 1; see DNA/SORT_NOTES.md).

## 7. File types

| Name | Meaning |
|---|---|
| `trajectory.dat / .dat.zst` | All saved frames (snapshots) of a simulation. Huge. Stored compressed as .zst. |
| `last_conf.dat / initial_conf.dat` | Single snapshot: the final / first configuration of a stage. |
| `energy.dat` | Energy values over time, one line per logged step. |
| `structure.top (.top)` | Topology: which nucleotide is bonded to which, and the sequence. |
| `input.txt` | Simulation settings (temperature, salt, steps, model). One per stage. |
| `simulation.log` | Text log written by the simulator. |
| `.zst` | zstd-compressed file. Needs `zstd -d file.zst` before use (see Read this before using). |
| `MANIFEST.csv` | List of every file in RNA/ or DNA/: path, size, date, design, multiplicity, run, role. |
| `SORT_NOTES.md` | Notes from the sorting step: where files came from, what is uncertain. |
| `README.md (in a run folder)` | Run summary: design, dates, parameters, status, which run replaces which. |
| `README_SORT.md (umbrella)` | Same idea for an umbrella-sampling source folder. |
| `ORIGINAL_NAMES.tsv` | Table: new file name and the original name it had on the source machine. |
| `window_NN/ (umbrella)` | One umbrella window: a simulation held at one value of the reaction coordinate, with equil/ and prod/ parts. |
| `WHAM profile (wham_profile.dat)` | Free-energy profile along the coordinate, combined from all windows by WHAM. |
| `parameters.txt / us_template_input.txt` | Umbrella set-up: window range, spring stiffness, step counts; template of the window input. |
| `.spm / .svg (AFM)` | .spm: raw AFM image file from the microscope. .svg: exported vector picture. |
| `run_NN` | Run number inside a design. run_01 is the newest main production run. |
| `monomer / 8mer / 3mer / design` | Number of units simulated together. 'design' = multiplicity not specified (nanotweezer). |

## 8. Full folder tree (depth up to 4, counts and sizes per folder)

<details>
<summary>RNA/  - 1962 files, 5.7 GB</summary>

- `RNA/` : 1962 files, 5.7 GB
  - `ffs/` : 544 files, 871.5 MB
    - `ffs_K1IVT8mer/` : 544 files, 871.5 MB
      - `21shoot` : 18 files, 85.1 MB
      - `32shoot` : 217 files, 130.4 MB
      - `43shoot` : 110 files, 60.1 MB
      - `54shoot` : 118 files, 59.2 MB
      - `HBKL5` : 72 files, 38.1 MB
  - `misc/` : 70 files, 80.4 MB
    - `SDoff_loose` : 3 files, 925.7 kB
    - `Sorting` : 3 files, 64.9 MB
    - `T3IVT_3mer/` : 14 files, 722.0 kB
      - `T3IVT_3mer_analysis` : 1 files, 9.1 kB
      - `T3IVT_3mer_md` : 1 files, 9.1 kB
      - `T3IVT_3mer_relax` : 9 files, 544.8 kB
    - `monomer_analysis/` : 10 files, 1.5 MB
      - `T3s` : 4 files, 809.6 kB
    - `plot_0527/` : 40 files, 12.3 MB
      - `K1IVT8mer` : 9 files, 901.3 kB
      - `T3IVT8mer` : 10 files, 8.9 MB
      - `mean_structures/` : 6 files, 881.6 kB
        - `K1IVT8mer` : 3 files, 427.7 kB
        - `T3IVT8mer` : 3 files, 453.9 kB
  - `simulations/` : 593 files, 4.1 GB
    - `2H_standard/` : 104 files, 1.2 GB
      - `8mer/` : 57 files, 983.4 MB
        - `run_01` : 57 files, 983.4 MB
      - `monomer/` : 47 files, 187.9 MB
        - `run_01` : 47 files, 187.9 MB
    - `K1_IVT_curved/` : 164 files, 1.1 GB
      - `8mer/` : 76 files, 901.2 MB
        - `run_01` : 62 files, 893.5 MB
        - `run_02` : 14 files, 7.6 MB
      - `monomer/` : 88 files, 238.4 MB
        - `run_01` : 53 files, 173.4 MB
        - `run_02` : 35 files, 65.0 MB
    - `K1_curved/` : 151 files, 1.2 GB
      - `8mer/` : 62 files, 913.2 MB
        - `run_01` : 46 files, 903.0 MB
        - `run_02` : 16 files, 10.2 MB
      - `monomer/` : 89 files, 316.5 MB
        - `run_01` : 51 files, 174.1 MB
        - `run_02` : 38 files, 142.4 MB
    - `T3_IVT_wavy/` : 87 files, 256.3 MB
      - `monomer/` : 87 files, 256.3 MB
        - `run_01` : 49 files, 188.9 MB
        - `run_02` : 38 files, 67.3 MB
    - `T3_wavy/` : 87 files, 255.5 MB
      - `monomer/` : 87 files, 255.5 MB
        - `run_01` : 51 files, 188.1 MB
        - `run_02` : 36 files, 67.4 MB
  - `starting_structures` : 14 files, 534.6 kB
  - `umbrella/` : 738 files, 728.1 MB
    - `umbrella_K1IVT8mer/` : 406 files, 44.8 MB
      - `umbrella/` : 392 files, 44.1 MB
        - `window_00` : 14 files, 1.7 MB
        - `window_01` : 14 files, 1.7 MB
        - `window_02` : 14 files, 1.7 MB
        - `window_03` : 14 files, 1.7 MB
        - `window_04` : 14 files, 1.7 MB
        - `window_05` : 14 files, 1.7 MB
        - `window_06` : 14 files, 1.7 MB
        - `window_07` : 14 files, 1.7 MB
        - `window_08` : 14 files, 1.7 MB
        - `window_09` : 14 files, 1.7 MB
        - `window_10` : 14 files, 1.7 MB
        - `window_11` : 14 files, 1.7 MB
        - `window_12` : 14 files, 1.7 MB
        - `window_13` : 14 files, 1.7 MB
        - `window_14` : 14 files, 1.7 MB
        - `window_15` : 14 files, 1.7 MB
        - `window_16` : 14 files, 1.7 MB
        - `window_17` : 14 files, 1.7 MB
        - `window_18` : 14 files, 1.7 MB
        - `window_19` : 14 files, 1.7 MB
        - `window_20` : 14 files, 1.7 MB
        - `window_21` : 14 files, 1.7 MB
        - `window_22` : 14 files, 1.7 MB
        - `window_23` : 14 files, 1.7 MB
        - `window_24` : 14 files, 1.7 MB
      - `umbrella_utils` : 9 files, 65.8 kB
    - `us_k1/` : 332 files, 683.3 MB
      - `umbrella/` : 317 files, 679.5 MB
        - `window_00` : 15 files, 33.9 MB
        - `window_01` : 15 files, 33.9 MB
        - `window_02` : 15 files, 33.9 MB
        - `window_03` : 15 files, 33.9 MB
        - `window_04` : 15 files, 33.9 MB
        - `window_05` : 15 files, 33.9 MB
        - `window_06` : 15 files, 33.9 MB
        - `window_07` : 15 files, 33.9 MB
        - `window_08` : 15 files, 33.9 MB
        - `window_09` : 15 files, 34.0 MB
        - `window_10` : 15 files, 33.9 MB
        - `window_11` : 15 files, 34.0 MB
        - `window_12` : 15 files, 33.9 MB
        - `window_13` : 15 files, 33.9 MB
        - `window_14` : 15 files, 34.0 MB
        - `window_15` : 15 files, 33.9 MB
        - `window_16` : 15 files, 34.0 MB
        - `window_17` : 15 files, 33.9 MB
        - `window_18` : 15 files, 33.9 MB
        - `window_19` : 15 files, 34.0 MB
      - `usutils` : 8 files, 65.3 kB

</details>

<details>
<summary>DNA/  - 3826 files, 12.3 GB</summary>

- `DNA/` : 3826 files, 12.3 GB
  - `.pytest_cache/` : 6 files, 955 B
    - `v/` : 3 files, 425 B
      - `cache` : 3 files, 425 B
  - `__pycache__` : 2 files, 103.0 kB
  - `nanotw_figures/` : 33 files, 15.1 MB
    - `angle_distribution` : 8 files, 370.0 kB
    - `fig_simulation/` : 22 files, 14.7 MB
      - `avg_close` : 3 files, 1.2 MB
      - `avg_open` : 3 files, 754.1 kB
      - `colored functional parts` : 11 files, 7.5 MB
  - `simulations/` : 262 files, 10.0 GB
    - `T3_IVT_wavy/` : 147 files, 6.2 GB
      - `3mer/` : 26 files, 63.1 MB
        - `run_01` : 26 files, 63.1 MB
      - `8mer/` : 70 files, 5.5 GB
        - `run_01` : 70 files, 5.5 GB
      - `monomer/` : 51 files, 693.2 MB
        - `run_01` : 51 files, 693.2 MB
    - `chol_fluoro_nanotw/` : 43 files, 196.1 MB
      - `design/` : 43 files, 196.1 MB
        - `run_01` : 23 files, 129.0 MB
        - `run_02` : 20 files, 67.2 MB
    - `nanotw_closed_patch4/` : 33 files, 439.8 MB
      - `design/` : 33 files, 439.8 MB
        - `run_01` : 33 files, 439.8 MB
    - `nanotw_open_patch4/` : 34 files, 3.1 GB
      - `design/` : 34 files, 3.1 GB
        - `run_01` : 34 files, 3.1 GB
    - `nanotw_starting_structures` : 5 files, 3.4 MB
  - `tests_sort/` : 4 files, 56.8 kB
    - `__pycache__` : 2 files, 47.6 kB
  - `umbrella/` : 3515 files, 2.3 GB
    - `nanotw_umbrella_sampling/` : 15 files, 7.9 MB
      - `prep` : 9 files, 4.2 MB
    - `nanotw_ustest/` : 415 files, 279.4 MB
      - `prep` : 9 files, 4.2 MB
      - `umbrella/` : 393 files, 267.8 MB
        - `db` : 86 files, 53.2 MB
        - `window_00` : 15 files, 21.9 MB
        - `window_01` : 11 files, 7.5 MB
        - `window_02` : 11 files, 7.5 MB
        - `window_03` : 11 files, 7.5 MB
        - `window_04` : 11 files, 7.5 MB
        - `window_05` : 11 files, 7.5 MB
        - `window_06` : 11 files, 7.5 MB
        - `window_07` : 11 files, 7.5 MB
        - `window_08` : 11 files, 7.5 MB
        - `window_09` : 11 files, 7.5 MB
        - `window_10` : 12 files, 11.1 MB
        - `window_11` : 11 files, 7.5 MB
        - `window_12` : 11 files, 7.5 MB
        - `window_13` : 11 files, 7.5 MB
        - `window_14` : 11 files, 7.5 MB
        - `window_15` : 11 files, 7.6 MB
        - `window_16` : 11 files, 7.6 MB
        - `window_17` : 11 files, 7.6 MB
        - `window_18` : 11 files, 7.5 MB
        - `window_19` : 11 files, 7.5 MB
        - `window_20` : 11 files, 7.6 MB
        - `window_21` : 11 files, 7.5 MB
        - `window_22` : 11 files, 7.5 MB
        - `window_23` : 11 files, 7.5 MB
        - `window_24` : 33 files, 15.4 MB
    - `umbrella_DNA/` : 1339 files, 698.1 MB
      - `umbrella/` : 509 files, 342.8 MB
        - `db` : 86 files, 53.2 MB
        - `window_00` : 14 files, 11.1 MB
        - `window_01` : 14 files, 11.1 MB
        - `window_02` : 14 files, 11.1 MB
        - `window_03` : 14 files, 11.1 MB
        - `window_04` : 14 files, 11.1 MB
        - `window_05` : 14 files, 11.1 MB
        - `window_06` : 14 files, 11.1 MB
        - `window_07` : 14 files, 11.1 MB
        - `window_08` : 14 files, 11.1 MB
        - `window_09` : 14 files, 11.1 MB
        - `window_10` : 14 files, 11.1 MB
        - `window_11` : 14 files, 11.1 MB
        - `window_12` : 14 files, 11.1 MB
        - `window_13` : 14 files, 11.1 MB
        - `window_14` : 14 files, 11.1 MB
        - `window_15` : 14 files, 11.1 MB
        - `window_16` : 14 files, 11.1 MB
        - `window_17` : 14 files, 11.1 MB
        - `window_18` : 14 files, 11.1 MB
        - `window_19` : 14 files, 11.1 MB
        - `window_20` : 14 files, 11.1 MB
        - `window_21` : 14 files, 11.1 MB
        - `window_22` : 14 files, 11.1 MB
        - `window_23` : 14 files, 11.1 MB
        - `window_24` : 36 files, 18.9 MB
      - `umbrella_500confs/` : 540 files, 343.8 MB
        - `db` : 86 files, 53.2 MB
        - `window_00` : 14 files, 11.2 MB
        - `window_01` : 14 files, 11.2 MB
        - `window_02` : 14 files, 11.2 MB
        - `window_03` : 14 files, 11.2 MB
        - `window_04` : 14 files, 11.2 MB
        - `window_05` : 14 files, 11.2 MB
        - `window_06` : 14 files, 11.2 MB
        - `window_07` : 14 files, 11.2 MB
        - `window_08` : 14 files, 11.2 MB
        - `window_09` : 14 files, 11.2 MB
        - `window_10` : 14 files, 11.2 MB
        - `window_11` : 14 files, 11.2 MB
        - `window_12` : 14 files, 11.2 MB
        - `window_13` : 14 files, 11.2 MB
        - `window_14` : 14 files, 11.2 MB
        - `window_15` : 14 files, 11.2 MB
        - `window_16` : 14 files, 11.2 MB
        - `window_17` : 14 files, 11.2 MB
        - `window_18` : 14 files, 11.2 MB
        - `window_19` : 14 files, 11.2 MB
        - `window_20` : 14 files, 11.2 MB
        - `window_21` : 14 files, 11.2 MB
        - `window_22` : 14 files, 11.2 MB
        - `window_23` : 14 files, 11.2 MB
        - `window_24` : 36 files, 19.0 MB
      - `umbrella_k1/` : 129 files, 5.2 MB
        - `umbrella` : 122 files, 1.6 MB
      - `us_k1/` : 158 files, 6.3 MB
        - `umbrella` : 151 files, 2.7 MB
    - `umbrella_DNA_opening/` : 435 files, 284.9 MB
      - `umbrella/` : 420 files, 281.1 MB
        - `window_00` : 14 files, 11.2 MB
        - `window_01` : 14 files, 11.2 MB
        - `window_02` : 14 files, 11.2 MB
        - `window_03` : 14 files, 11.2 MB
        - `window_04` : 14 files, 11.2 MB
        - `window_05` : 14 files, 11.2 MB
        - `window_06` : 14 files, 11.2 MB
        - `window_07` : 14 files, 11.2 MB
        - `window_08` : 14 files, 11.2 MB
        - `window_09` : 14 files, 11.2 MB
        - `window_10` : 14 files, 11.2 MB
        - `window_11` : 14 files, 11.2 MB
        - `window_12` : 14 files, 11.2 MB
        - `window_13` : 14 files, 11.2 MB
        - `window_14` : 14 files, 11.2 MB
        - `window_15` : 14 files, 11.2 MB
        - `window_16` : 14 files, 11.2 MB
        - `window_17` : 14 files, 11.2 MB
        - `window_18` : 14 files, 11.2 MB
        - `window_19` : 14 files, 11.2 MB
        - `window_20` : 14 files, 11.2 MB
        - `window_21` : 14 files, 11.2 MB
        - `window_22` : 14 files, 11.2 MB
        - `window_23` : 14 files, 11.2 MB
        - `window_24` : 14 files, 11.2 MB
      - `utils` : 3 files, 16.1 kB
    - `umbrella_K1IVT8mer/` : 411 files, 45.0 MB
      - `umbrella/` : 396 files, 44.4 MB
        - `window_00` : 14 files, 1.7 MB
        - `window_01` : 14 files, 1.7 MB
        - `window_02` : 14 files, 1.7 MB
        - `window_03` : 14 files, 1.7 MB
        - `window_04` : 14 files, 1.7 MB
        - `window_05` : 14 files, 1.7 MB
        - `window_06` : 14 files, 1.7 MB
        - `window_07` : 14 files, 1.7 MB
        - `window_08` : 14 files, 1.7 MB
        - `window_09` : 14 files, 1.7 MB
        - `window_10` : 14 files, 1.7 MB
        - `window_11` : 14 files, 1.7 MB
        - `window_12` : 14 files, 1.7 MB
        - `window_13` : 14 files, 1.7 MB
        - `window_14` : 14 files, 1.7 MB
        - `window_15` : 14 files, 1.7 MB
        - `window_16` : 14 files, 1.7 MB
        - `window_17` : 14 files, 1.7 MB
        - `window_18` : 14 files, 1.7 MB
        - `window_19` : 14 files, 1.7 MB
        - `window_20` : 14 files, 1.7 MB
        - `window_21` : 14 files, 1.7 MB
        - `window_22` : 14 files, 1.7 MB
        - `window_23` : 14 files, 1.7 MB
        - `window_24` : 14 files, 1.7 MB
      - `umbrella_utils` : 9 files, 65.8 kB
    - `umbrella_T3IVT/` : 467 files, 897.6 MB
      - `umbrella/` : 459 files, 897.5 MB
        - `window_00` : 15 files, 35.9 MB
        - `window_01` : 15 files, 35.9 MB
        - `window_02` : 15 files, 35.9 MB
        - `window_03` : 15 files, 36.0 MB
        - `window_04` : 15 files, 35.9 MB
        - `window_05` : 15 files, 35.8 MB
        - `window_06` : 15 files, 35.9 MB
        - `window_07` : 15 files, 35.8 MB
        - `window_08` : 15 files, 35.8 MB
        - `window_09` : 15 files, 35.9 MB
        - `window_10` : 15 files, 35.9 MB
        - `window_11` : 15 files, 35.8 MB
        - `window_12` : 15 files, 35.8 MB
        - `window_13` : 15 files, 35.8 MB
        - `window_14` : 15 files, 35.8 MB
        - `window_15` : 15 files, 35.8 MB
        - `window_16` : 15 files, 35.9 MB
        - `window_17` : 15 files, 35.8 MB
        - `window_18` : 15 files, 35.8 MB
        - `window_19` : 15 files, 35.8 MB
        - `window_20` : 15 files, 35.8 MB
        - `window_21` : 15 files, 35.8 MB
        - `window_22` : 15 files, 35.8 MB
        - `window_23` : 15 files, 35.8 MB
        - `window_24` : 15 files, 35.8 MB
      - `usutil` : 2 files, 1.0 kB
    - `umbrella_T3IVT8mer/` : 433 files, 47.1 MB
      - `umbrella/` : 427 files, 46.5 MB
        - `window_00` : 14 files, 1.8 MB
        - `window_01` : 14 files, 1.8 MB
        - `window_02` : 14 files, 1.8 MB
        - `window_03` : 14 files, 1.8 MB
        - `window_04` : 14 files, 1.8 MB
        - `window_05` : 14 files, 1.8 MB
        - `window_06` : 14 files, 1.8 MB
        - `window_07` : 14 files, 1.8 MB
        - `window_08` : 14 files, 1.8 MB
        - `window_09` : 14 files, 1.8 MB
        - `window_10` : 14 files, 1.8 MB
        - `window_11` : 14 files, 1.8 MB
        - `window_12` : 14 files, 1.8 MB
        - `window_13` : 14 files, 1.8 MB
        - `window_14` : 14 files, 1.8 MB
        - `window_15` : 14 files, 1.8 MB
        - `window_16` : 14 files, 1.8 MB
        - `window_17` : 14 files, 1.8 MB
        - `window_18` : 14 files, 1.8 MB
        - `window_19` : 14 files, 1.8 MB
        - `window_20` : 14 files, 1.8 MB
        - `window_21` : 14 files, 1.8 MB
        - `window_22` : 14 files, 1.8 MB
        - `window_23` : 14 files, 1.8 MB
        - `window_24` : 14 files, 1.8 MB

</details>

<details>
<summary>WetLab/  - 98 files, 347.0 MB</summary>

- `WetLab/` : 98 files, 347.0 MB
  - `AFM/` : 98 files, 347.0 MB
    - `20260826` : 1 files, 850.0 kB
    - `20260828` : 10 files, 14.9 MB
    - `20260901` : 3 files, 7.9 MB
    - `20260902` : 15 files, 28.5 MB
    - `20260903` : 27 files, 95.8 MB
    - `20260907` : 12 files, 26.0 MB
    - `20260908` : 9 files, 65.4 MB
    - `20260909` : 3 files, 7.3 MB
    - `20260910` : 18 files, 100.2 MB

</details>

Questions about scripts: repo docs/RNA.md and repo docs/DNA.md.

