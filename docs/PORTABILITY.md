# Portability: hard-coded paths and cluster settings

No path was edited (scientific logic untouched). To run elsewhere, replace these prefixes (`sed -i`) in the scripts you use.

Occurrences of absolute-path prefixes (count, prefix):

```
  98 /home/dfan/oxDNA
  80 /home/dfan/rna_ys
  75 /Users/a000000/db
  50 /home/donglifan/oxDNA
  10 /home/dfan/inputs_paras_oxDNA
   4 /home/donglifan/oxdna
   3 /home/dfan/rna_utils
   1 /home/donglifan/rna_utils
   1 /home/donglifan/inputs_paras_oxDNA
   1 /home/donglifan/dna_utils
   1 /Users/yunshiyang/Downloads
```

Main locations: `rna/umbrella/run_all.sh`, `dna/umbrella/run_all.sh`, `run_production_only.sh`, `run_prod.sh`, `rna/thermal_ramp/run_thermal_ramp.sh`, `setup_run.py`, `setup_analysis.py`, `setup_kl_ffs.py`, `setup_prod.py`, `thermal_ramp_setup.py`, `ffs_*.py`, `setup_equil_cmds.py`, and `templates/`.

- `/home/dfan/oxDNA/build/bin/oxDNA` and `/home/donglifan/oxDNA/...`: the oxDNA binary (see README, "Installing oxDNA").
- `/home/dfan/rna_utils`, `/home/dfan/inputs_paras_oxDNA`: template directories (now in `templates/`).
- `/Users/a000000/db/...` (dna/umbrella/run_all.sh): local mount of the cluster; replace.
- `/Users/yunshiyang/Downloads` (notebooks/rna/figure_template.ipynb): local output folder.

## Slurm (jobs/)
- Partitions: `gpu`, `gpu-a100-small`, `compute`; account `research-as-bn`; these are site-specific.
- Jobs run `conda activate oxdna` (an env that holds oxDNA's Python bindings; differs from `yunshi-oxdna` in environment.yml).
- Array pattern: `sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cmds` runs line N of a `cmds` file per task.

## DelftBlue (setup_hpc.sh)

Facts below come from the job scripts in `jobs/` and read-only checks on the cluster (`module spider`, `sinfo`).

| Item | Value | Source |
|---|---|---|
| Account | `research-as-bn` | all `jobs/*/*.sh` |
| Partitions used | `gpu-a100-small` (production arrays, 4 h), `gpu` (MD, 24 h), `compute` (MC / CPU, 1 h) | `jobs/` |
| GPU request | `--gpus-per-task=1` | `jobs/` |
| Env in jobs | `conda activate oxdna` (miniconda3 in `/home/donglifan`) | `jobs/` |
| oxDNA binary on DelftBlue | `/home/donglifan/oxdna/build/bin/oxDNA` (lowercase `oxdna`) | `dna/setup/setup_run.py`, `setup_equil_cmds.py` |
| `/home/dfan/oxDNA/...` | path from the other (WSL) machine, not DelftBlue | most RNA scripts |
| Modules used by `setup_hpc.sh` | `2025` stack, then `gcc/12.4.0 cmake/3.30.5 cuda/12.5 python/3.11.9` | exist per `module spider`; load checked on a login node |

The existing job scripts load no modules (they rely on the conda env). `setup_hpc.sh` uses a venv instead, so its jobs need `module load ...` and `source $ENV_ROOT/bin/activate`; no new Slurm templates were added. Settings at the top of `setup_hpc.sh`: `SCRATCH_BASE`, `ENV_ROOT`, `OXDNA_DIR`, `MODULES`.

Tested: `module load` of the module set above (login node, read-only). Not tested: running `setup_hpc.sh`, the CUDA build, `cmake`/`make` of oxDNA, and the local `--build-oxdna` option.
