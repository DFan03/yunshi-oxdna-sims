# yunshi-oxdna-sims

Analysis, setup and job scripts for coarse-grained oxDNA/oxRNA simulations: RNA kissing-loop designs (2H, K1, T3, IVT; monomer and 8mer) and DNA nanotweezers (nanotw, T3IVT) with umbrella sampling. Scripts only; simulation data are not in git.

## Quick start

```bash
git clone https://github.com/DFan03/yunshi-oxdna-sims.git
cd yunshi-oxdna-sims
bash setup.sh        # conda/mamba env if available, else .venv + pip; runs a smoke test
```

`setup.sh` installs numpy, scipy, pandas, matplotlib, seaborn, jupyter and `oxDNA_analysis_tools` (oat), then imports them and byte-compiles every script. It warns (not fails) if `oxpy` is missing.

## Guides

- [docs/RNA.md](docs/RNA.md): RNA workflow, script by script (when to use, in/out, example command, versions).
- [docs/DNA.md](docs/DNA.md): same for DNA.
- [docs/SCRIPTS.md](docs/SCRIPTS.md): short index of scripts by stage.
- [docs/VERSIONS.md](docs/VERSIONS.md): canonical vs archived versions of each script.

Also: [docs/PORTABILITY.md](docs/PORTABILITY.md) (hard-coded paths, Slurm settings). Self-contained HTML twins `docs/RNA.html` and `docs/DNA.html` exist for local browsing (GitHub shows them as source).

## Repo map

| Folder | Contents |
|---|---|
| `rna/` | `setup/`, `analysis/`, `multimer/`, `ffs/`, `thermal_ramp/`, `umbrella/` |
| `dna/` | `setup/`, `analysis/`, `multimer/`, `thermal_ramp/`, `umbrella/` |
| `templates/` | oxDNA input templates (pre-relax, CUDA relax, MD, umbrella, FFS) |
| `jobs/` | Slurm scripts |
| `notebooks/` | plotting notebooks (outputs stripped) |
| `run_configs/` | per-run configs that differ from the templates |
| `archive/` | older distinct versions of scripts |
| `tools/` | helper tools (notebook stripper, docs generator) |

## How the scripts chain

`setup_run.py` -> MC pre-relax (CPU) -> CUDA MD relax -> production MD -> analysis (oat, `hb_cal.py`, bend/twist) and/or umbrella sampling + WHAM, FFS (RNA only), thermal ramp -> notebooks. Run scripts inside the simulation/analysis directory.

## Data location

Not in git. Data live on the USB drive in `Yunshi_data/` (`RNA/`, `DNA/`, `WetLab/`); see `FILES.md` there for the inventory.

## Versions policy

Identical copies were merged by md5; the newest modification time is canonical and sits in `rna/` or `dna/`. Older distinct versions are in `archive/<name>/<yyyymmdd>_<tag>.<ext>`. `docs/VERSIONS.md` says which is canonical and what differs. Per-window `check_equil.py` files (auto-generated, differ only in ANCHOR_A/ANCHOR_B/XI0) are summarised in `archive/check_equil/variants_*.tsv`.

## Installing oxDNA (external prerequisite)

The simulator itself is not included. Scripts call `/home/dfan/oxDNA/build/bin/oxDNA` and need the Python bindings (`oxpy`, not on PyPI) for `hb_cal`, `helix_bt`, `multimer_tb`, `ffs_*`, `persistence_length`. Version used: commit `8eed395` (3 Mar 2026).

```bash
git clone https://github.com/lorenzo-rovigatti/oxDNA.git && cd oxDNA
git checkout 8eed395357533d0a81aecfb73a0f228aab4f180c
mkdir build && cd build
cmake .. -DPython=ON            # add -DCUDA=ON for GPU (needs CUDA toolkit)
make -j$(nproc) && make install # oxpy installed into the active Python env
```

Then edit the binary path (see `docs/PORTABILITY.md`). Build options: see the oxDNA documentation.

## Notes

Scripts were copied unmodified; hard-coded cluster paths remain (documented). Notebook outputs are stripped by `tools/strip_notebook.py`.
