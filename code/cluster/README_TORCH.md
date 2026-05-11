# Torch Cluster Guide

This folder contains the active Torch launcher for the Python direct-geometry
calibration.

Source of truth for live calibration state:

```bash
CALIBRATION_STATUS.md
```

Active model code:

```bash
code/model/
```

## Account Rule

Always use:

```bash
torch_pr_570_general
```

Do not use `torch_pr_571_general` for this project.

## Where To Submit

On Torch, submit jobs from:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/cluster
```

The active submission script is:

```bash
submit_python_direct_geometry_overnight.sh
```

The script uses relative paths to `../model`, so run it from this folder.

## Smoke Test

Use a small array and short evaluation budget before any overnight run:

```bash
DT_DIRECT_RUN_TAG=py_direct_smoke_YYYYMMDD \
DT_DIRECT_SETUP=fast \
DT_DIRECT_MAX_ITER_EQ=4 \
DT_DIRECT_MAX_EVALS=2 \
DT_DIRECT_BUDGET_SEC=1800 \
sbatch --array=1-2%2 submit_python_direct_geometry_overnight.sh
```

## Overnight Run

The short partition accepts jobs up to `03:55:00`. A standard wide run is:

```bash
DT_DIRECT_RUN_TAG=py_direct_renewal_calibrated_global_YYYYMMDD \
DT_DIRECT_SETUP=benchmark \
DT_DIRECT_BOUNDS=global \
sbatch --array=1-40%40 submit_python_direct_geometry_overnight.sh
```

## Monitor

```bash
squeue -u td2248
tail -f logs/slurm_py_direct_JOB_TASK.out
```

The active script writes worker outputs under:

```bash
results_python_direct_geometry_${DT_DIRECT_RUN_TAG}/
```

## Collect Results

Use the Python collector from the port:

```bash
cd $SCRATCH/projects/Fertility_Spring26/code/model
python tools/collect_direct_geometry_results.py \
  --results-dir ../cluster/results_python_direct_geometry_RUN_TAG
```

Historical MATLAB launchers and old result directories are archived in
`calibration_archive/cluster_matlab_2026-05-07/`.
