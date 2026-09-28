#!/usr/bin/env bash
#SBATCH --job-name=e5f_price_diag_tests
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/fertility_night_calibration_20260928_v1/price_diagnostic_v1
sha256sum -c source_manifest.sha256
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python -m unittest -v test_diagnose_e5f_evening_housing_failures
