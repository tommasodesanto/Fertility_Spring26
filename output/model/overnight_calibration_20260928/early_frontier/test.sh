#!/usr/bin/env bash
#SBATCH --job-name=e5f_frontier_tests
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
cd "$stage/early_frontier_v1"
sha256sum -c source_manifest.sha256
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=1
export PYTHONPATH="$PWD:$original/code/model/tools"
export REVIEWED_RECOVERY_SEARCH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4/run_e5f_utility_comparison_search.py"
apptainer exec --bind "$stage/project:$original" --pwd "$PWD" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python -m unittest -v test_e5f_early_fertility_frontier
