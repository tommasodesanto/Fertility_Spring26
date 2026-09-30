#!/usr/bin/env bash
#SBATCH --job-name=e5f_fixed_credit_runtime_check
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
set -euo pipefail
[[ "${1:-}" == "--self-test-only" ]] || { echo 'prepared only: production requires lead review and staged source pin' >&2; exit 2; }
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/fixed_credit_runtime_validation_20260929
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
source=$remote/source_runtime_validation_v1
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=/work/numba_cache
mkdir -p "$remote/results" "$remote/cache"
[[ -f "$source/plan.json" && -f "$source/run_runtime_validation.py" ]] || { echo 'staged source absent' >&2; exit 2; }
[[ ! -e "$remote/results/controller_mock_tests_v1" ]] || { echo 'refusing duplicate output' >&2; exit 2; }
apptainer exec --bind "$reference:$original:ro,$source:/work/fixed_credit_source:ro,$remote/results:/work/results,$remote/cache:/work/numba_cache" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/fixed_credit_source/run_runtime_validation.py --plan /work/fixed_credit_source/plan.json --output /work/results/controller_mock_tests_v1 --self-test-only
