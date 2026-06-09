# Cluster Calibration

Active launchers:

- `submit_python_direct_geometry_overnight.sh`
- `submit_intergen_housing_fertility_twohour_panel.sh`
- `submit_intergen_housing_fertility_global_de.sh`

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

Two-hour one-market global-DE panel:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/cluster
INTERGEN_RUN_TAG=intergen_candidate_no_timing_v0_globalde_3g_20260609 \
INTERGEN_GLOBAL_EVALS_PER_TASK=320 \
INTERGEN_GLOBAL_POP_SIZE=24 \
sbatch --array=1-64%64 --time=2:05:00 --mem=3G \
  submit_intergen_housing_fertility_global_de.sh
```

This runs `intergen_housing_fertility.cli global-de-panel`: independent
Latin-hypercube restarts followed by differential-evolution proposals over the
same 13 economic parameters used by `candidate_no_timing_v0`. It uses the same
model solution, grid, income-state count, target set, and checkpoint format as
the two-hour panel launcher; only the parameter proposal algorithm changes.

After the tasks finish:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/model
python tools/collect_intergen_panel_results.py \
  --results-dir ../cluster/results_intergen_housing_fertility_intergen_candidate_no_timing_v0_twohour_20260608
```

Torch shutdown snapshot workflow, June 2026:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/cluster
./collect_intergen_shutdown_snapshot.sh snapshot_rome0650_20260609
```

The snapshot script collects all matching one-market result directories,
writes `combined_top200.csv` and `combined_summary.json`, and archives the
snapshot under `shutdown_snapshots/`. The optional local pull helper copies
those tarballs back to `output/model/intergen_shutdown_snapshots/`:

```bash
code/cluster/pull_intergen_shutdown_snapshots_local.sh 2026-06-09T06:42:00 25 45
```

Active collected result:

- `results_python_direct_geometry_py_direct_renewal_calibrated_global_12h_20260506/`

Historical MATLAB launchers and older result directories were archived on
2026-05-07 under `calibration_archive/cluster_matlab_2026-05-07/`.
