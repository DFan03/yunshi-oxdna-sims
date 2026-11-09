# Script versions

How the scripts were collected: every script-like file (`*.py *.sh *.ipynb`, job scripts, `cmds*`, `input*`) found in the two machine trees (`wsl:/home/dfan` = RNA work, `db:/scratch/donglifan` = DNA work, plus the local staged copies) was hashed (md5; notebooks are hashed **after output stripping**). Identical files were collapsed. Where several different versions of one script name exist **within the same family (rna / dna)**, the version with the **newest modification time** is canonical and lives in the main tree; every other distinct version is kept under `archive/<name>/<yyyymmdd>_<family>-<machine>-<source-dir>.<ext>`. Nothing was dropped silently.

Rules and exceptions:

- **Canonical = newest mtime** among distinct versions of the same name in the same family. mtime is the only evidence available for "which run used it"; where a name has several versions the table lists the directories (runs) that used each one.
- `check_equil.py` is auto-generated per umbrella window by `setup_umb.py` (it differs only by 3 constants `ANCHOR_A`, `ANCHOR_B`, `XI0`). One canonical copy is kept; **all variants are listed** in `archive/check_equil/variants_rna.tsv` and `variants_dna.tsv` (source dir, constants, md5) instead of storing hundreds of near-identical files.
- Per-run configuration files (`input_pre_relax`, `input_cuda_relax`, `input_md`, `cmds*`, `parameters.txt`, `forces.txt`, `us_template_input.txt`, `op.txt`, `helixces.txt`, `kissingloop.csv`, ...) are *run settings, not code versions*. The generic template (the one in `rna_utils/` or `dna_utils/`) is canonical in `templates/<family>/`; every run-specific file that differs from it is kept at its original relative path under `run_configs/<family>/` (tree mirrors the original machine tree). Files identical to the template are not repeated; see counts below.
- The same script content often appears in both families (e.g. `hb_cal.py` in `rna/analysis/` and `dna/analysis/`; `dna_utils` was copied from `rna_utils`). They are kept in both trees so each family is self-contained; the md5 is identical when listed as "same as rna/..." below.
- Notebooks have outputs and execution counts removed (`tools/strip_notebook.py`), so different notebook "versions" that only differed by cell outputs are collapsed.
- The two machines have different oxDNA paths (`/home/dfan/oxDNA` vs `/home/donglifan/oxdna`), which is why some DNA/RNA copies of the same script are distinct versions (see `docs/PORTABILITY.md`).

## Canonical scripts (main tree)

| script | family | canonical file | mtime | newest source dir | copies (dirs) | distinct versions | archived? |
|---|---|---|---|---|---|---|---|
| `angle_cal.py` | dna | `dna/analysis/angle_cal.py` | 20260609 | `nanotw/close_path4_20c_seqD/close_path4_20c_seqD_analysis` | 1 (1) | 1 |  |
| `calculate_bond.py` | dna | `dna/analysis/calculate_bond.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `calculate_bond_instability.py` | dna | `dna/analysis/calculate_bond_instability.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `dis_bent_plotter.py` | dna | `dna/analysis/dis_bent_plotter.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `dis_bond_plotter.py` | dna | `dna/analysis/dis_bond_plotter.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `dis_energy_plotter.py` | dna | `dna/analysis/dis_energy_plotter.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `dis_twist_plotter.py` | dna | `dna/analysis/dis_twist_plotter.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `distance_cal.py` | dna | `dna/analysis/distance_cal.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `distance_two_sets.py` | dna | `dna/analysis/distance_two_sets.py` | 20260414 | `nanotw/ustest/prep` | 2 (2) | 1 |  |
| `hb_cal.py` | dna | `dna/analysis/hb_cal.py` | 20260617 | `dna_utils/analysis_tools` | 8 (4) | 1 |  |
| `helical_bent_twist_cal.py` | dna | `dna/analysis/helical_bent_twist_cal.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `helix_bt.py` | dna | `dna/analysis/helix_bt.py` | 20260617 | `dna_utils/analysis_tools` | 4 (2) | 1 |  |
| `plot_energydat.py` | dna | `dna/analysis/plot_energydat.py` | 20260617 | `dna_utils` | 2 (1) | 1 |  |
| `signed_helical_bent_twist_cal.py` | dna | `dna/analysis/signed_helical_bent_twist_cal.py` | 20260413 | `T3IVT/T3IVT_analysis` | 2 (1) | 1 |  |
| `tb2.py` | dna | `dna/analysis/tb2.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `total_hbond_cal.py` | dna | `dna/analysis/total_hbond_cal.py` | 20260617 | `dna_utils/analysis_tools` | 4 (2) | 1 |  |
| `twist_bent.py` | dna | `dna/analysis/twist_bent.py` | 20260617 | `dna_utils/analysis_tools` | 2 (1) | 1 |  |
| `multimer_tb.py` | dna | `dna/multimer/multimer_tb.py` | 20260508 | `T3IVT_8mer/T3IVT_8mer_md` | 1 (1) | 1 |  |
| `multimer_forces.py` | dna | `dna/setup/multimer_forces.py` | 20260506 | `T3IVT_8mer` | 2 (1) | 1 |  |
| `prepare_oat_workflow.py` | dna | `dna/setup/prepare_oat_workflow.py` | 20260617 | `dna_utils` | 2 (1) | 1 |  |
| `replacer.py` | dna | `dna/setup/replacer.py` | 20260617 | `dna_utils` | 2 (1) | 1 |  |
| `setup_analysis.py` | dna | `dna/setup/setup_analysis.py` | 20260511 | `T3IVT/T3IVT_analysis` | 2 (1) | 1 |  |
| `setup_equil_cmds.py` | dna | `dna/setup/setup_equil_cmds.py` | 20260715 | `umbrella_DNA_opening/umbrella` | 11 (6) | 2 | yes |
| `setup_run.py` | dna | `dna/setup/setup_run.py` | 20260617 | `chol_fluoro_nanotw` | 11 (8) | 3 | yes |
| `thermal_ramp_setup.py` | dna | `dna/thermal_ramp/thermal_ramp_setup.py` | 20260617 | `dna_utils/thermal_ramp` | 2 (1) | 1 |  |
| `wrap_thermal_ramp.py` | dna | `dna/thermal_ramp/wrap_thermal_ramp.py` | 20260617 | `dna_utils/thermal_ramp` | 2 (1) | 1 |  |
| `2state.py` | dna | `dna/umbrella/2state.py` | 20260617 | `umbrella_T3IVT8mer/umbrella` | 4 (2) | 1 |  |
| `check_dis.py` | dna | `dna/umbrella/check_dis.py` | 20260707 | `umbrella_DNA_opening` | 3 (2) | 1 |  |
| `check_equil.py` | dna | `dna/umbrella/check_equil.py` | 20260715 | `nanotw/ustest/umbrella/db/umbrella/window_04` | 302 (180) | 100 | yes |
| `compute_barrier_forces.py` | dna | `dna/umbrella/compute_barrier_forces.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 10 (5) | 1 |  |
| `delta_energy.py` | dna | `dna/umbrella/delta_energy.py` | 20260617 | `umbrella_T3IVT/umbrella` | 2 (1) | 1 |  |
| `force_from_pmf.py` | dna | `dna/umbrella/force_from_pmf.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 10 (5) | 1 |  |
| `plot_com_by_window.py` | dna | `dna/umbrella/plot_com_by_window.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 6 (3) | 1 |  |
| `plot_complete.py` | dna | `dna/umbrella/plot_complete.py` | 20260609 | `umbrella_T3IVT8mer/umbrella` | 4 (2) | 1 |  |
| `plot_completed.py` | dna | `dna/umbrella/plot_completed.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 12 (7) | 2 | yes |
| `plot_window.py` | dna | `dna/umbrella/plot_window.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 6 (3) | 1 |  |
| `pwindow.py` | dna | `dna/umbrella/pwindow.py` | 20260609 | `umbrella_T3IVT8mer/umbrella` | 4 (2) | 1 |  |
| `run_all.sh` | dna | `dna/umbrella/run_all.sh` | 20260715 | `umbrella_DNA_opening/utils` | 22 (13) | 2 | yes |
| `run_production_only.sh` | dna | `dna/umbrella/run_production_only.sh` | 20260422 | `umbrella_DNA/umbrella_500confs` | 8 (5) | 2 | yes |
| `setup_prod.py` | dna | `dna/umbrella/setup_prod.py` | 20260715 | `umbrella_DNA_opening/umbrella` | 6 (3) | 1 |  |
| `setup_umb.py` | dna | `dna/umbrella/setup_umb.py` | 20260707 | `umbrella_DNA_opening` | 7 (4) | 2 | yes |
| `setup_umbrella.py` | dna | `dna/umbrella/setup_umbrella.py` | 20260414 | `nanotw/umbrella_sampling` | 1 (1) | 1 |  |
| `state_average.py` | dna | `dna/umbrella/state_average.py` | 20260618 | `umbrella_T3IVT8mer/umbrella` | 4 (2) | 1 |  |
| `tb_wham.py` | dna | `dna/umbrella/tb_wham.py` | 20260618 | `umbrella_T3IVT/umbrella` | 2 (1) | 1 |  |
| `try.sh` | dna | `dna/umbrella/try.sh` | 20260422 | `umbrella_DNA/umbrella_500confs` | 4 (2) | 1 |  |
| `tune_stiff.py` | dna | `dna/umbrella/tune_stiff.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 26 (15) | 2 | yes |
| `wham_cal.py` | dna | `dna/umbrella/wham_cal.py` | 20260609 | `umbrella_T3IVT8mer/umbrella` | 4 (2) | 1 |  |
| `wham_calculate.py` | dna | `dna/umbrella/wham_calculate.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 6 (3) | 1 |  |
| `wham_com.py` | dna | `dna/umbrella/wham_com.py` | 20260724 | `umbrella_DNA_opening/umbrella` | 6 (3) | 1 |  |
| `cpu_mono.sh` | dna | `jobs/dna/cpu_mono.sh` | 20260617 | `umbrella_T3IVT/umbrella` | 2 (1) | 1 |  |
| `equil_job.sh` | dna | `jobs/dna/equil_job.sh` | 20260715 | `umbrella_DNA_opening/umbrella` | 9 (5) | 1 |  |
| `mcjob.sh` | dna | `jobs/dna/mcjob.sh` | 20260507 | `T3IVT/T3IVT_relax` | 5 (4) | 1 |  |
| `mdjob.sh` | dna | `jobs/dna/mdjob.sh` | 20260507 | `T3IVT/T3IVT_md` | 8 (6) | 2 | yes |
| `prod_job.sh` | dna | `jobs/dna/prod_job.sh` | 20260715 | `umbrella_DNA_opening/umbrella` | 10 (5) | 3 | yes |
| `queue_25_prod.sh` | dna | `jobs/dna/queue_25_prod.sh` | 20260422 | `umbrella_DNA/umbrella_500confs` | 2 (1) | 1 |  |
| `queue_relax.sh` | dna | `jobs/dna/queue_relax.sh` | 20260617 | `chol_fluoro_nanotw/chol_fluoro_nanotw_relax` | 1 (1) | 1 |  |
| `queue.sh` | rna | `jobs/rna/queue.sh` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella` | 2 (1) | 1 |  |
| `queue_25_prod.sh` | rna | `jobs/rna/queue_25_prod.sh` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 5 (4) | 1 |  |
| `calculate_bond.py` | rna | `rna/analysis/calculate_bond.py` | 20260318 | `tests/SD_MR_tests/2H_control/2H_control_37c_seqdOff_analysis` | 9 (9) | 1 |  |
| `calculate_bond_instability.py` | rna | `rna/analysis/calculate_bond_instability.py` | 20260318 | `ests/SD_MR_tests/K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_relax` | 20 (20) | 1 |  |
| `dis_bent_plotter.py` | rna | `rna/analysis/dis_bent_plotter.py` | 20260323 | `s/SD_MR_tests/K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_analysis` | 7 (7) | 3 | yes |
| `dis_bond_plotter.py` | rna | `rna/analysis/dis_bond_plotter.py` | 20260323 | `ts/SD_MR_tests/K1IVT_remodeled/K1IVT_remodeled_SM20_analysis` | 8 (8) | 4 | yes |
| `dis_energy_plotter.py` | rna | `rna/analysis/dis_energy_plotter.py` | 20260318 | `tests/SD_MR_tests/2H_control/2H_control_37c_seqdOff_analysis` | 8 (8) | 2 | yes |
| `dis_twist_plotter.py` | rna | `rna/analysis/dis_twist_plotter.py` | 20260323 | `s/SD_MR_tests/K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_analysis` | 7 (7) | 3 | yes |
| `distance_cal.py` | rna | `rna/analysis/distance_cal.py` | 20260318 | `tests/SD_MR_tests/2H_control/2H_control_37c_seqdOff_analysis` | 8 (8) | 1 |  |
| `hb_cal.py` | rna | `rna/analysis/hb_cal.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 21 (19) | 1 |  |
| `helical_bent_twist_cal.py` | rna | `rna/analysis/helical_bent_twist_cal.py` | 20260318 | `tests/SD_MR_tests/2H_control/2H_control_37c_seqdOff_analysis` | 5 (5) | 1 |  |
| `helix_bt.py` | rna | `rna/analysis/helix_bt.py` | 20260318 | `rna_utils/analysis_tools` | 13 (12) | 1 |  |
| `log_transformer.py` | rna | `rna/analysis/log_transformer.py` | 20260311 | `rna_ys/default_runs/T3_IVT/oxdna_simulated` | 1 (1) | 1 |  |
| `persistence_length.py` | rna | `rna/analysis/persistence_length.py` | 20260521 | `SDoff_MRRon_RUN/K1_20c_sdoff_mron/K1_20c_sdoff_mron_analysis` | 2 (2) | 2 | yes |
| `plot_bond_instability.py` | rna | `rna/analysis/plot_bond_instability.py` | 20260318 | `ests/SD_MR_tests/K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_relax` | 8 (8) | 1 |  |
| `plot_energydat.py` | rna | `rna/analysis/plot_energydat.py` | 20260315 | `relax_21md_seqd_mmr/K1IVT_37relax_21md_seqd_mmr_thermal_ramp` | 3 (2) | 1 |  |
| `signed_helical_bent_twist_cal.py` | rna | `rna/analysis/signed_helical_bent_twist_cal.py` | 20260330 | `rna_utils/analysis_tools` | 7 (6) | 1 |  |
| `tb2.py` | rna | `rna/analysis/tb2.py` | 20260318 | `tests/SD_MR_tests/2H_control/2H_control_37c_seqdOff_analysis` | 6 (6) | 1 |  |
| `total_hbond_cal.py` | rna | `rna/analysis/total_hbond_cal.py` | 20260318 | `rna_utils/analysis_tools` | 12 (11) | 1 |  |
| `twist_bent.py` | rna | `rna/analysis/twist_bent.py` | 20260318 | `tests/SD_MR_tests/2H_control/2H_control_37c_seqdOff_analysis` | 6 (6) | 1 |  |
| `extract_conf_hb.py` | rna | `rna/ffs/extract_conf_hb.py` | 20260407 | `rna_ys/SDoff_MRRon_RUN/ffs_K1IVT8mer` | 2 (2) | 1 |  |
| `ffs_flux.py` | rna | `rna/ffs/ffs_flux.py` | 20260325 | `rna_ys/tests/ffs_final_test/FLUX` | 5 (3) | 1 |  |
| `ffs_shoot.py` | rna | `rna/ffs/ffs_shoot.py` | 20260408 | `rna_ys/SDoff_MRRon_RUN/ffs_K1IVT8mer/21shoot` | 15 (14) | 10 | yes |
| `setup_kl_ffs.py` | rna | `rna/ffs/setup_kl_ffs.py` | 20260325 | `rna_ys/tests/ffs_final_test` | 1 (1) | 1 |  |
| `multimer_tb.py` | rna | `rna/multimer/multimer_tb.py` | 20260331 | `off_MRRon_RUN/2H8mer_20c_sdoff_mron/2H8mer_20c_sdoff_mron_md` | 14 (13) | 3 | yes |
| `pairwise_tb.py` | rna | `rna/multimer/pairwise_tb.py` | 20260812 | `Ron_RUN/K1IVT8mer_20c_sdoff_mron/K1IVT8mer_20c_sdoff_mron_md` | 2 (2) | 1 |  |
| `al_mt_creater.py` | rna | `rna/setup/al_mt_creater.py` | 20260326 | `rna_ys/tests/SD_MR_tests/alternative8merk1ivt` | 1 (1) | 1 |  |
| `force_remover.py` | rna | `rna/setup/force_remover.py` | 20260324 | `a_ys/tests/SD_MR_tests/K1IVT_remodeled/K1IVT_remodeled_relax` | 1 (1) | 1 |  |
| `force_remover.py.txt` | rna | `rna/setup/force_remover.py.txt` | 20260324 | `rna_utils` | 1 (1) | 1 |  |
| `mt_creater.py` | rna | `rna/setup/mt_creater.py` | 20260326 | `rna_ys/tests/SD_MR_tests/alternative8merk1ivt` | 1 (1) | 1 |  |
| `multimerforce.py` | rna | `rna/setup/multimerforce.py` | 20260329 | `_MRRon_RUN/K18mer_20c_sdoff_mron/K18mer_20c_sdoff_mron_relax` | 1 (1) | 1 |  |
| `ntremover.py` | rna | `rna/setup/ntremover.py` | 20260316 | `ys/tests/K1IVT_mutual_trap_test/K1IVT_mutual_trap_test_relax` | 1 (1) | 1 |  |
| `prepare_oat_workflow.py` | rna | `rna/setup/prepare_oat_workflow.py` | 20260318 | `ests/SD_MR_tests/K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_relax` | 20 (19) | 2 | yes |
| `replacer.py` | rna | `rna/setup/replacer.py` | 20260318 | `ests/SD_MR_tests/K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_relax` | 11 (10) | 1 |  |
| `setup_analysis.py` | rna | `rna/setup/setup_analysis.py` | 20260421 | `MRRon_RUN/T3IVT_20c_sdoff_mron/T3IVT_20c_sdoff_mron_analysis` | 11 (10) | 2 | yes |
| `setup_run.py` | rna | `rna/setup/setup_run.py` | 20260423 | `rna_ys/SDoff_MRRon_RUN/T3IVT_3mer` | 28 (27) | 2 | yes |
| `special_replacer.py` | rna | `rna/setup/special_replacer.py` | 20260310 | `rna_ys/default_runs/K1IVT_8mer` | 2 (2) | 1 |  |
| `run_thermal_ramp.sh` | rna | `rna/thermal_ramp/run_thermal_ramp.sh` | 20260321 | `s/SD_MR_tests/K1IVT_8mer_mtrap/K1IVT_8mer_mtrap_thermal_ramp` | 2 (2) | 2 | yes |
| `thermal_ramp_setup.py` | rna | `rna/thermal_ramp/thermal_ramp_setup.py` | 20260321 | `rna_ys/tests/SD_MR_tests/K1IVT_8mer_mtrap` | 4 (3) | 2 | yes |
| `wrap_thermal_ramp.py` | rna | `rna/thermal_ramp/wrap_thermal_ramp.py` | 20260323 | `s/SD_MR_tests/K1IVT_8mer_mtrap/K1IVT_8mer_mtrap_thermal_ramp` | 4 (3) | 1 |  |
| `2state.py` | rna | `rna/umbrella/2state.py` | 20260617 | `umbrella_K1IVT8mer/umbrella` | 2 (1) | 1 |  |
| `check_dis.py` | rna | `rna/umbrella/check_dis.py` | 20260422 | `umbrella_DNA/umbrella_k1` | 5 (3) | 1 |  |
| `check_equil.py` | rna | `rna/umbrella/check_equil.py` | 20260507 | `rna_ys/umbrella_K1IVT8mer/umbrella/window_00` | 127 (90) | 53 | yes |
| `compute_barrier_forces.py` | rna | `rna/umbrella/compute_barrier_forces.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 9 (7) | 1 |  |
| `force_from_pmf.py` | rna | `rna/umbrella/force_from_pmf.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 9 (7) | 1 |  |
| `plot_com_by_window.py` | rna | `rna/umbrella/plot_com_by_window.py` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella` | 2 (1) | 1 |  |
| `plot_complete.py` | rna | `rna/umbrella/plot_complete.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 7 (6) | 1 |  |
| `plot_completed.py` | rna | `rna/umbrella/plot_completed.py` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella` | 2 (1) | 1 |  |
| `plot_window.py` | rna | `rna/umbrella/plot_window.py` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella` | 2 (1) | 1 |  |
| `pwindow.py` | rna | `rna/umbrella/pwindow.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 7 (6) | 1 |  |
| `run_all.sh` | rna | `rna/umbrella/run_all.sh` | 20260422 | `rna_ys/us_k1/umbrella` | 5 (3) | 2 | yes |
| `run_prod.sh` | rna | `rna/umbrella/run_prod.sh` | 20260422 | `rna_ys/us_k1/umbrella` | 1 (1) | 1 |  |
| `run_production_only.sh` | rna | `rna/umbrella/run_production_only.sh` | 20260422 | `rna_ys/us_k1/umbrella` | 3 (2) | 2 | yes |
| `setup_prod.py` | rna | `rna/umbrella/setup_prod.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 7 (6) | 1 |  |
| `setup_umb.py` | rna | `rna/umbrella/setup_umb.py` | 20260507 | `umbrella_K1IVT8mer` | 8 (5) | 2 | yes |
| `state_average.py` | rna | `rna/umbrella/state_average.py` | 20260618 | `umbrella_K1IVT8mer/umbrella` | 2 (1) | 1 |  |
| `tune_stiff.py` | rna | `rna/umbrella/tune_stiff.py` | 20260511 | `umbrella_K1IVT8mer/umbrella` | 9 (7) | 2 | yes |
| `wham_cal.py` | rna | `rna/umbrella/wham_cal.py` | 20260520 | `umbrella_K1IVT8mer/umbrella` | 7 (6) | 2 | yes |
| `wham_calculate.py` | rna | `rna/umbrella/wham_calculate.py` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella` | 2 (1) | 1 |  |
| `wham_com.py` | rna | `rna/umbrella/wham_com.py` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella` | 2 (1) | 1 |  |
| `Untitled.ipynb` | dna | `notebooks/dna/Untitled.ipynb` | 20260415 | `nanotw/ustest` | 1 (1) | 1 |  |
| `plot_analysis_results.ipynb` | dna | `notebooks/dna/plot_analysis_results.ipynb` | 20260511 | `T3IVT/T3IVT_analysis` | 2 (1) | 1 |  |
| `plot_curves.ipynb` | dna | `notebooks/dna/plot_curves.ipynb` | 20260609 | `nanotw/angle_distribution` | 1 (1) | 1 |  |
| `plot_energy.ipynb` | dna | `notebooks/dna/plot_energy.ipynb` | 20260617 | `chol_fluoro_nanotw/chol_fluoro_nanotw_relax` | 10 (8) | 2 | yes |
| `plot_multimer.ipynb` | dna | `notebooks/dna/plot_multimer.ipynb` | 20260508 | `T3IVT_8mer/T3IVT_8mer_md` | 1 (1) | 1 |  |
| `Untitled.ipynb` | rna | `notebooks/rna/Untitled.ipynb` | 20260812 | `off_MRRon_RUN/2H8mer_20c_sdoff_mron/2H8mer_20c_sdoff_mron_md` | 8 (8) | 4 | yes |
| `cross_plotter.ipynb` | rna | `notebooks/rna/cross_plotter.ipynb` | 20260422 | `rna_ys/SDoff_MRRon_RUN/monomer_analysis/T3s` | 2 (2) | 1 |  |
| `figure_template.ipynb` | rna | `notebooks/rna/figure_template.ipynb` | 20260527 | `MRRon_RUN/K1IVT_20c_sdoff_mron/K1IVT_20c_sdoff_mron_analysis` | 1 (1) | 1 |  |
| `mutimer_plotter-abs.ipynb` | rna | `notebooks/rna/mutimer_plotter-abs.ipynb` | 20260527 | `rna_ys/plot_0527` | 1 (1) | 1 |  |
| `mutimer_plotter.ipynb` | rna | `notebooks/rna/mutimer_plotter.ipynb` | 20260527 | `rna_ys/plot_0527` | 1 (1) | 1 |  |
| `plot_abs.ipynb` | rna | `notebooks/rna/plot_abs.ipynb` | 20260527 | `rna_ys/plot_0527` | 1 (1) | 1 |  |
| `plot_analysis_results.ipynb` | rna | `notebooks/rna/plot_analysis_results.ipynb` | 20260527 | `MRRon_RUN/K1IVT_20c_sdoff_mron/K1IVT_20c_sdoff_mron_analysis` | 9 (8) | 4 | yes |
| `plot_energy.ipynb` | rna | `notebooks/rna/plot_energy.ipynb` | 20260423 | `rna_ys/SDoff_MRRon_RUN/T3IVT_3mer/T3IVT_3mer_relax` | 39 (38) | 3 | yes |
| `plot_monos.ipynb` | rna | `notebooks/rna/plot_monos.ipynb` | 20260527 | `rna_ys/plot_0527` | 1 (1) | 1 |  |
| `plot_multimer.ipynb` | rna | `notebooks/rna/plot_multimer.ipynb` | 20260420 | `Ron_RUN/K1IVT8mer_20c_sdoff_mron/K1IVT8mer_20c_sdoff_mron_md` | 7 (6) | 5 | yes |

## Archived older / alternative versions

Short diff = line counts and first changed lines of the archived version relative to the canonical one (`+` lines are in the archived file, `-` lines are in the canonical file).

| script | family | archived file | mtime | used in (example dir) | #paths | what differs |
|---|---|---|---|---|---|---|
| `Untitled.ipynb` | rna | `archive/Untitled/20260310_rna-wsl-rna_ys-tests-analysis_test.ipynb` | 20260310 | `rna_ys/tests/analysis_test/Untitled.ipynb` | 1 | notebook source differs |
| `Untitled.ipynb` | rna | `archive/Untitled/20260331_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_8mer_mtrap-K1IVT_8mer_mtrap_md37.ipynb` | 20260331 | `_tests/K1IVT_8mer_mtrap/K1IVT_8mer_mtrap_md37/Untitled.ipynb` | 5 | notebook source differs |
| `Untitled.ipynb` | rna | `archive/Untitled/20260812_rna-wsl-s-SDoff_MRRon_RUN-K1IVT8mer_20c_sdoff_mron-K1IVT8mer_20c_sdoff_mron_md.ipynb` | 20260812 | `er_20c_sdoff_mron/K1IVT8mer_20c_sdoff_mron_md/Untitled.ipynb` | 1 | notebook source differs |
| `dis_bent_plotter.py` | rna | `archive/dis_bent_plotter/20260318_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_mmr_seqd_37c-K1IVT_mmr_seqd_37c_md.py` | 20260318 | `K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_md/dis_bent_plotter.py` | 1 | +4/-4 lines vs canonical; e.g. `-    fig, ax1 = plt.subplots(figsize=(15, 5))` | `+    fig, ax1 = plt.subplots(figsize=(30, 5))` |
| `dis_bent_plotter.py` | rna | `archive/dis_bent_plotter/20260318_rna-wsl-rna_ys-tests-distance_analysis_test_K1IVT-utils.py` | 20260318 | `tests/distance_analysis_test_K1IVT/utils/dis_bent_plotter.py` | 2 | +4/-6 lines vs canonical; e.g. `-OUTPUT_FILE = "kl_distance_vs_helical_bent.png"  ` | `-    fig, ax1 = plt.subplots(figsize=(15, 5))` |
| `dis_bond_plotter.py` | rna | `archive/dis_bond_plotter/20260318_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_mmr_seqd_37c-K1IVT_mmr_seqd_37c_md.py` | 20260318 | `K1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_md/dis_bond_plotter.py` | 1 | +7/-8 lines vs canonical; e.g. `-OUTPUT_FILE = "kl_distance_vs_KL_HBond.png" ` | `+OUTPUT_FILE = "kl_distance_vs_total_HBond.png" ` |
| `dis_bond_plotter.py` | rna | `archive/dis_bond_plotter/20260318_rna-wsl-rna_ys-tests-distance_analysis_test_K1IVT-utils.py` | 20260318 | `tests/distance_analysis_test_K1IVT/utils/dis_bond_plotter.py` | 3 | +7/-9 lines vs canonical; e.g. `-OUTPUT_FILE = "kl_distance_vs_KL_HBond.png" ` | `+` |
| `dis_bond_plotter.py` | rna | `archive/dis_bond_plotter/20260323_rna-wsl-rna_ys-tests-SD_MR_tests-2H_control-2H_control_37c_seqdOff_analysis.py` | 20260323 | `_control/2H_control_37c_seqdOff_analysis/dis_bond_plotter.py` | 1 | +1/-1 lines vs canonical; e.g. `-    ax2.set_ylim([-0.5,7])` | `+    ax2.set_ylim([0,7])` |
| `dis_energy_plotter.py` | rna | `archive/dis_energy_plotter/20260316_rna-wsl-rna_ys-tests-distance_analysis_test_K1IVT-37c_K1_IVT.py` | 20260316 | `istance_analysis_test_K1IVT/37c_K1_IVT/dis_energy_plotter.py` | 1 | +1/-1 lines vs canonical; e.g. `-            energy = float(parts[1])` | `+            energy = float(parts[3])` |
| `dis_twist_plotter.py` | rna | `archive/dis_twist_plotter/20260318_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_mmr_seqd_37c-K1IVT_mmr_seqd_37c_md.py` | 20260318 | `1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_md/dis_twist_plotter.py` | 1 | +4/-4 lines vs canonical; e.g. `-    fig, ax1 = plt.subplots(figsize=(15, 5))` | `+    fig, ax1 = plt.subplots(figsize=(30, 5))` |
| `dis_twist_plotter.py` | rna | `archive/dis_twist_plotter/20260318_rna-wsl-rna_ys-tests-distance_analysis_test_K1IVT-utils.py` | 20260318 | `ests/distance_analysis_test_K1IVT/utils/dis_twist_plotter.py` | 2 | +5/-6 lines vs canonical; e.g. `-OUTPUT_FILE = "kl_distance_vs_helical_twist.png" ` | `+` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260325_rna-wsl-rna_ys-tests-ffs_final_test-10_shoot.py` | 20260325 | `rna_ys/tests/ffs_final_test/10_shoot/ffs_shoot.py` | 1 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../21_shoot/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260325_rna-wsl-rna_ys-tests-ffs_final_test-32_shoot.py` | 20260325 | `rna_ys/tests/ffs_final_test/32_shoot/ffs_shoot.py` | 2 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../43_shoot/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260325_rna-wsl-rna_ys-tests-ffs_final_test-43_shoot.py` | 20260325 | `rna_ys/tests/ffs_final_test/43_shoot/ffs_shoot.py` | 2 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../54_shoot/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260325_rna-wsl-rna_ys-tests-ffs_final_test-54_shoot.py` | 20260325 | `rna_ys/tests/ffs_final_test/54_shoot/ffs_shoot.py` | 3 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../FLUX/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260325_rna-wsl-rna_ys-tests-ffs_from_trj_test-21_shoot.py` | 20260325 | `rna_ys/tests/ffs_from_trj_test/21_shoot/ffs_shoot.py` | 2 | +1/-1 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../32_shoot/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260325_rna-wsl-rna_ys-tests-ffs_from_trj_test-54_shoot.py` | 20260325 | `rna_ys/tests/ffs_from_trj_test/54_shoot/ffs_shoot.py` | 1 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../HBKL6/*.dat'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260407_rna-wsl-rna_ys-SDoff_MRRon_RUN-ffs_K1IVT8mer-32shoot.py` | 20260407 | `rna_ys/SDoff_MRRon_RUN/ffs_K1IVT8mer/32shoot/ffs_shoot.py` | 1 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../43shoot/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260407_rna-wsl-rna_ys-SDoff_MRRon_RUN-ffs_K1IVT8mer-43shoot.py` | 20260407 | `rna_ys/SDoff_MRRon_RUN/ffs_K1IVT8mer/43shoot/ffs_shoot.py` | 1 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../54shoot/success*'` |
| `ffs_shoot.py` | rna | `archive/ffs_shoot/20260407_rna-wsl-rna_ys-SDoff_MRRon_RUN-ffs_K1IVT8mer-54shoot.py` | 20260407 | `rna_ys/SDoff_MRRon_RUN/ffs_K1IVT8mer/54shoot/ffs_shoot.py` | 1 | +4/-4 lines vs canonical; e.g. `-starting_conf_pattern = '../32shoot/success*'` | `+starting_conf_pattern = '../HBKL5/conf*'` |
| `mdjob.sh` | dna | `archive/mdjob/20260507_dna-db-T3IVT-T3IVT_relax.sh` | 20260507 | `T3IVT/T3IVT_relax/mdjob.sh` | 5 | +2/-2 lines vs canonical; e.g. `-#SBATCH --partition=gpu` | `-#SBATCH --time=24:00:00` |
| `multimer_tb.py` | rna | `archive/multimer_tb/20260331_rna-wsl-rna_ys-SDoff_MRRon_RUN-K1IVT8mer_20c_sdoff_mron.py` | 20260331 | `a_ys/SDoff_MRRon_RUN/K1IVT8mer_20c_sdoff_mron/multimer_tb.py` | 5 | +91/-237 lines vs canonical; e.g. `-` | `+import numpy as np` |
| `multimer_tb.py` | rna | `archive/multimer_tb/20260331_rna-wsl-rna_ys-SDoff_MRRon_RUN.py` | 20260331 | `rna_ys/SDoff_MRRon_RUN/multimer_tb.py` | 4 | +62/-204 lines vs canonical; e.g. `-` | `+import numpy as np` |
| `persistence_length.py` | rna | `archive/persistence_length/20260521_rna-wsl-_ys-SDoff_MRRon_RUN-K1IVT_20c_sdoff_mron-K1IVT_20c_sdoff_mron_analysis.py` | 20260521 | `off_mron/K1IVT_20c_sdoff_mron_analysis/persistence_length.py` | 1 | +34/-61 lines vs canonical; e.g. `-persistence_length_csv.py - PROPER STANDARD PRESENTATION` | `-` |
| `plot_analysis_results.ipynb` | rna | `archive/plot_analysis_results/20260330_rna-wsl-rna_ys-SDoff_MRRon_RUN-2H_20c_sdoff_mron-2H_20c_sdoff_mron_analysis.ipynb` | 20260330 | `_mron/2H_20c_sdoff_mron_analysis/plot_analysis_results.ipynb` | 1 | notebook source differs |
| `plot_analysis_results.ipynb` | rna | `archive/plot_analysis_results/20260422_rna-wsl-_ys-SDoff_MRRon_RUN-T3IVT_20c_sdoff_mron-T3IVT_20c_sdoff_mron_analysis.ipynb` | 20260422 | `on/T3IVT_20c_sdoff_mron_analysis/plot_analysis_results.ipynb` | 2 | notebook source differs |
| `plot_analysis_results.ipynb` | rna | `archive/plot_analysis_results/20260422_rna-wsl-rna_ys-SDoff_MRRon_RUN-T3_20c_sdoff_mron-T3_20c_sdoff_mron_analysis.ipynb` | 20260422 | `_mron/T3_20c_sdoff_mron_analysis/plot_analysis_results.ipynb` | 5 | notebook source differs |
| `plot_completed.py` | dna | `archive/plot_completed/20260422_dna-db-umbrella_DNA-umbrella_500confs-window_24.py` | 20260422 | `umbrella_DNA/umbrella_500confs/window_24/plot_completed.py` | 6 | +11/-11 lines vs canonical; e.g. `-    Violin plot of COM–COM distance per completed window.` | `-    """` |
| `plot_energy.ipynb` | rna | `archive/plot_energy/20260309_rna-wsl-rna_ys-default_runs-K1_IVT.ipynb` | 20260309 | `rna_ys/default_runs/K1_IVT/plot_energy.ipynb` | 1 | notebook source differs |
| `plot_energy.ipynb` | rna | `archive/plot_energy/20260318_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_mmr_seqd_37c-K1IVT_mmr_seqd_37c_relax.ipynb` | 20260318 | `1IVT_mmr_seqd_37c/K1IVT_mmr_seqd_37c_relax/plot_energy.ipynb` | 10 | notebook source differs |
| `plot_energy.ipynb` | dna | `archive/plot_energy/20260320_dna-wsl-nanotw-open_patch4_20c_seqD-open_patch4_20c_seqD_relax.ipynb` | 20260320 | `patch4_20c_seqD/open_patch4_20c_seqD_relax/plot_energy.ipynb` | 1 | notebook source differs |
| `plot_multimer.ipynb` | rna | `archive/plot_multimer/20260331_rna-wsl-rna_ys-SDoff_MRRon_RUN-2H8mer_20c_sdoff_mron-2H8mer_20c_sdoff_mron_md.ipynb` | 20260331 | `_20c_sdoff_mron/2H8mer_20c_sdoff_mron_md/plot_multimer.ipynb` | 3 | notebook source differs |
| `plot_multimer.ipynb` | rna | `archive/plot_multimer/20260331_rna-wsl-rna_ys-SDoff_MRRon_RUN-K18mer_20c_sdoff_mron-K18mer_20c_sdoff_mron_md.ipynb` | 20260331 | `_20c_sdoff_mron/K18mer_20c_sdoff_mron_md/plot_multimer.ipynb` | 1 | notebook source differs |
| `plot_multimer.ipynb` | rna | `archive/plot_multimer/20260331_rna-wsl-rna_ys-SDoff_MRRon_RUN.ipynb` | 20260331 | `rna_ys/SDoff_MRRon_RUN/plot_multimer.ipynb` | 1 | notebook source differs |
| `plot_multimer.ipynb` | rna | `archive/plot_multimer/20260331_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_8mer_mtrap-K1IVT_8mer_mtrap_md37.ipynb` | 20260331 | `s/K1IVT_8mer_mtrap/K1IVT_8mer_mtrap_md37/plot_multimer.ipynb` | 1 | notebook source differs |
| `prepare_oat_workflow.py` | rna | `archive/prepare_oat_workflow/20260318_rna-wsl-sts-K1_IVT_LOWMATCH_SeqD_MMR_37c-K1_IVT_LOWMATCH_SeqD_MMR_37c_analysis.py` | 20260318 | `1_IVT_LOWMATCH_SeqD_MMR_37c_analysis/prepare_oat_workflow.py` | 8 | +1/-1 lines vs canonical; e.g. `-src_conf = "../last_conf.dat"` | `+src_conf = "../mc_relax_last_conf.dat"` |
| `prod_job.sh` | dna | `archive/prod_job/20260422_dna-db-umbrella_DNA-umbrella_500confs.sh` | 20260422 | `umbrella_DNA/umbrella_500confs/prod_job.sh` | 4 | +4/-13 lines vs canonical; e.g. `+#` | `-#SBATCH --time=4:00:00` |
| `prod_job.sh` | dna | `archive/prod_job/20260609_dna-db-umbrella_T3IVT8mer-umbrella.sh` | 20260609 | `umbrella_T3IVT8mer/umbrella/prod_job.sh` | 4 | +1/-1 lines vs canonical; e.g. `-CMD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cmds)` | `+CMD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" prod_cmds)` |
| `run_all.sh` | dna | `archive/run_all/20260422_dna-db-umbrella_DNA-umbrella_500confs-db-umbrella.sh` | 20260422 | `umbrella_DNA/umbrella_500confs/db/umbrella/run_all.sh` | 20 | +174/-174 lines vs canonical; e.g. `-# window_00  xi0 = 9.0000` | `-echo '[window_00] equil starting (xi0=9.000)...'` |
| `run_all.sh` | rna | `archive/run_all/20260422_rna-db-umbrella_DNA-umbrella_k1-umbrella.sh` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella/run_all.sh` | 4 | +27/-137 lines vs canonical; e.g. `-# window_01  xi0 = 7.6421` | `-echo '[window_01] equil starting (xi0=7.642)...'` |
| `run_production_only.sh` | dna | `archive/run_production_only/20260422_dna-db-umbrella_DNA-umbrella_500confs-window_24.sh` | 20260422 | `rella_DNA/umbrella_500confs/window_24/run_production_only.sh` | 6 | +25/-25 lines vs canonical; e.g. `-cd window_00/prod && oxDNA input.txt > run.log 2>&1` | `+cd window_00/prod && /home/dfan/oxDNA/build/bin/oxDNA input` |
| `run_production_only.sh` | rna | `archive/run_production_only/20260422_rna-db-umbrella_DNA-umbrella_k1-umbrella.sh` | 20260422 | `umbrella_DNA/umbrella_k1/umbrella/run_production_only.sh` | 2 | +10/-50 lines vs canonical; e.g. `-cd window_00/prod && /home/dfan/oxDNA/build/bin/oxDNA input` | `+cd window_00/prod && /home/donglifan/oxdna/build/bin/oxDNA ` |
| `run_thermal_ramp.sh` | rna | `archive/run_thermal_ramp/20260315_rna-wsl-s-K1IVT_37relax_21md_seqd_mmr-K1IVT_37relax_21md_seqd_mmr_thermal_ramp.sh` | 20260315 | `K1IVT_37relax_21md_seqd_mmr_thermal_ramp/run_thermal_ramp.sh` | 1 | +15/-9 lines vs canonical; e.g. `-/home/dfan/oxDNA/build/bin/oxDNA inputMD_step2_33c` | `+/home/dfan/oxDNA/build/bin/oxDNA inputMD_step2_36c` |
| `setup_analysis.py` | rna | `archive/setup_analysis/20260326_rna-wsl-rna_ys-tests-SD_MR_tests-K1IVT_remodeled-K1IVT_remodeled_SM20_analysis.py` | 20260326 | `VT_remodeled/K1IVT_remodeled_SM20_analysis/setup_analysis.py` | 4 | +7/-8 lines vs canonical; e.g. `-sleep 1` | `-` |
| `setup_equil_cmds.py` | dna | `archive/setup_equil_cmds/20260609_dna-db-umbrella_T3IVT8mer-umbrella.py` | 20260609 | `umbrella_T3IVT8mer/umbrella/setup_equil_cmds.py` | 2 | +1/-1 lines vs canonical; e.g. `-with open("cmds", "w") as f:` | `+with open("equil_cmds", "w") as f:` |
| `setup_run.py` | rna | `archive/setup_run/20260318_rna-wsl-rna_ys-tests-SD_MR_tests-K1_mmr_seqd_37c-K1_mmr_seqd_37c_relax.py` | 20260318 | `_MR_tests/K1_mmr_seqd_37c/K1_mmr_seqd_37c_relax/setup_run.py` | 10 | +2/-0 lines vs canonical; e.g. `+    shutil.copy(os.path.join(utils_dir, "prepare_oat_workfl` | `+    shutil.copy(os.path.join(utils_dir, "calculate_bond_ins` |
| `setup_run.py` | dna | `archive/setup_run/20260506_dna-db-T3IVT_8mer.py` | 20260506 | `T3IVT_8mer/setup_run.py` | 5 | +1/-1 lines vs canonical; e.g. `-    utils_dir = "/home/donglifan/dna_utils"` | `+    utils_dir = "/home/donglifan/rna_utils"` |
| `setup_run.py` | dna | `archive/setup_run/20260617_dna-wsl-nanotw-chol_fluoro_nanotw.py` | 20260617 | `nanotw/chol_fluoro_nanotw/setup_run.py` | 4 | +3/-3 lines vs canonical; e.g. `-    utils_dir = "/home/donglifan/dna_utils"` | `+    utils_dir = "/home/dfan/dna_utils"` |
| `setup_umb.py` | rna | `archive/setup_umb/20260422_rna-db-umbrella_DNA-umbrella_k1.py` | 20260422 | `umbrella_DNA/umbrella_k1/setup_umb.py` | 5 | +30/-1 lines vs canonical; e.g. `+        equil_cmds.append(f"# {win_name}  xi0 = {xi0:.4f}")` | `+        equil_cmds.append(f"echo '[{win_name}] equil starti` |
| `setup_umb.py` | dna | `archive/setup_umb/20260609_dna-db-umbrella_T3IVT8mer.py` | 20260609 | `umbrella_T3IVT8mer/setup_umb.py` | 4 | +36/-97 lines vs canonical; e.g. `-Distance metric: both anchor-anchor distance AND COM-COM of` | `+Distance metric: anchor-anchor distance (also used as the r` |
| `thermal_ramp_setup.py` | rna | `archive/thermal_ramp_setup/20260314_rna-wsl-rna_ys-tests-thermal_ramps-K1IVT_37relax_21md_seqd_mmr.py` | 20260314 | `rmal_ramps/K1IVT_37relax_21md_seqd_mmr/thermal_ramp_setup.py` | 1 | +11/-0 lines vs canonical; e.g. `+        if i == 1:` | `+            # find original steps line` |
| `tune_stiff.py` | rna | `archive/tune_stiff/20260507_rna-db-umbrella_K1IVT8mer-umbrella_utils.py` | 20260507 | `umbrella_K1IVT8mer/umbrella_utils/tune_stiff.py` | 6 | +4/-1 lines vs canonical; e.g. `-        lines.append(f"cd {prod_dir} && {oxdna_bin} input.t` | `+        lines.append(f'echo "[{w}] production starting..."'` |
| `tune_stiff.py` | dna | `archive/tune_stiff/20260609_dna-db-umbrella_T3IVT8mer-umbrella.py` | 20260609 | `umbrella_T3IVT8mer/umbrella/tune_stiff.py` | 4 | +1/-4 lines vs canonical; e.g. `-        lines.append(f'echo "[{w}] production starting..."'` | `-        lines.append(f"cd {prod_dir} && {oxdna_bin} input.t` |
| `wham_cal.py` | rna | `archive/wham_cal/20260507_rna-db-umbrella_K1IVT8mer-umbrella_utils.py` | 20260507 | `umbrella_K1IVT8mer/umbrella_utils/wham_cal.py` | 3 | +155/-210 lines vs canonical; e.g. `-      * Reads prod/trajectory.dat -> anchor-anchor distance` | `-      * Computes per-frame HB count using OAT hb_list (same` |

## Canonical files identical in the rna and dna trees

32 files are byte-identical copies:

`dna/umbrella/2state.py`, `dna/analysis/calculate_bond.py`, `dna/analysis/calculate_bond_instability.py`, `dna/umbrella/check_dis.py`, `dna/umbrella/compute_barrier_forces.py`, `dna/analysis/dis_energy_plotter.py`, `dna/analysis/distance_cal.py`, `dna/umbrella/force_from_pmf.py`, `dna/analysis/hb_cal.py`, `dna/analysis/helical_bent_twist_cal.py`, `dna/analysis/helix_bt.py`, `dna/multimer/multimer_tb.py`, `dna/umbrella/plot_com_by_window.py`, `dna/umbrella/plot_complete.py`, `dna/umbrella/plot_completed.py`, `notebooks/dna/plot_energy.ipynb`, `dna/analysis/plot_energydat.py`, `dna/umbrella/plot_window.py`, `dna/umbrella/pwindow.py`, `jobs/dna/queue_25_prod.sh`, `dna/setup/replacer.py`, `dna/umbrella/setup_prod.py`, `dna/analysis/signed_helical_bent_twist_cal.py`, `dna/umbrella/state_average.py`, `dna/analysis/tb2.py`, `dna/thermal_ramp/thermal_ramp_setup.py`, `dna/analysis/total_hbond_cal.py`, `dna/analysis/twist_bent.py`, `dna/umbrella/wham_cal.py`, `dna/umbrella/wham_calculate.py`, `dna/umbrella/wham_com.py`, `dna/thermal_ramp/wrap_thermal_ramp.py`


## Auto-generated per-window files

- `check_equil.py` (dna): 302 copies, 100 distinct (constants differ per window); table in `archive/check_equil/variants_dna.tsv`.
- `check_equil.py` (rna): 127 copies, 53 distinct (constants differ per window); table in `archive/check_equil/variants_rna.tsv`.

## Per-run configuration files

| name | family | copies | distinct versions | dirs | identical to template (copies) |
|---|---|---|---|---|---|
| `input_pre_relax` | rna | 45 | 36 | 44 | 2 |
| `input_cuda_relax` | rna | 38 | 35 | 37 | 2 |
| `cmds_relax` | rna | 26 | 1 | 26 | no template |
| `input` | rna | 16 | 5 | 15 | no template |
| `cmds.sh` | rna | 15 | 4 | 15 | no template |
| `input_md` | rna | 15 | 10 | 14 | 5 |
| `input_cuda_relax` | dna | 11 | 9 | 9 | 2 |
| `input_pre_relax` | dna | 11 | 9 | 9 | 2 |
| `cmds_md` | dna | 8 | 3 | 6 | no template |
| `cmds_relax` | dna | 8 | 2 | 7 | no template |
| `cmds` | dna | 7 | 3 | 4 | no template |
| `cmds` | rna | 6 | 3 | 4 | no template |
| `cmds_mc` | dna | 5 | 1 | 4 | no template |
| `parameters.txt` | dna | 3 | 3 | 3 | no template |
| `us_template_input.txt` | dna | 3 | 2 | 3 | 1 |
| `parameters.txt` | rna | 3 | 2 | 3 | no template |
| `us_template_input.txt` | rna | 3 | 1 | 3 | 3 |
| `cmds.sh` | dna | 2 | 1 | 1 | no template |
| `helixces.txt` | dna | 2 | 1 | 2 | no template |
| `input_md` | dna | 2 | 1 | 1 | 2 |
| `inputMD_step1_37c` | rna | 2 | 2 | 2 | no template |
| `readme` | rna | 2 | 2 | 2 | no template |
| `forces.txt` | dna | 1 | 1 | 1 | no template |
| `indexes_connect` | dna | 1 | 1 | 1 | no template |
| `kissingloop.csv` | dna | 1 | 1 | 1 | no template |
| `README` | rna | 1 | 1 | 1 | no template |
| `Readme.txt` | rna | 1 | 1 | 1 | no template |
| `helixces.txt` | rna | 1 | 1 | 1 | no template |
| `inputMD_step10_28c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step10_4c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step11_27c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step12_26c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step13_25c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step2_33c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step2_36c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step3_30c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step3_35c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step4_26c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step4_34c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step5_22c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step5_33c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step6_19c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step6_32c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step7_15c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step7_31c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step8_11c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step8_30c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step9_29c` | rna | 1 | 1 | 1 | no template |
| `inputMD_step9_8c` | rna | 1 | 1 | 1 | no template |
| `kissingloop.csv` | rna | 1 | 1 | 1 | no template |
| `op.txt` | rna | 1 | 1 | 1 | no template |
| `readme.txt` | rna | 1 | 1 | 1 | no template |
