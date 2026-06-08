# Cluster Calibration

Active launchers:

- `submit_python_direct_geometry_overnight.sh`
- `submit_intergen_housing_fertility_twohour_panel.sh`

Two-hour one-market intergen panel test:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/cluster
INTERGEN_RUN_TAG=intergen_candidate_no_timing_v0_twohour_20260608 \
sbatch --array=1-8%8 submit_intergen_housing_fertility_twohour_panel.sh
```

This runs `intergen_housing_fertility.cli local-panel` with
`target_set=candidate_no_timing_v0`, `J=16`, `Nb=60`, five income states, six
owner rungs, one worker per task, and a 115-minute internal budget inside a
2:10 SLURM walltime. Task 1 includes the deterministic anchor cases; other
tasks use random draws only.

After the tasks finish:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/model
python tools/collect_intergen_panel_results.py \
  --results-dir ../cluster/results_intergen_housing_fertility_intergen_candidate_no_timing_v0_twohour_20260608
```

Active collected result:

- `results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/`

Historical MATLAB launchers and older result directories were archived on
2026-05-07 under `calibration_archive/cluster_matlab_2026-05-07/`.
