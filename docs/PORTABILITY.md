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
