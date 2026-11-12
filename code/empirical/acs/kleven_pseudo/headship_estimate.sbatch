#!/usr/bin/env bash
# Unchanged national estimator on one pooled headship design (2026-11-12).
#SBATCH --job-name=acs-headship-est
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_headship_20261112/logs/est_%x_%j.log
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DESIGN="${DESIGN:?}" HEADSHIP_OUT="${HEADSHIP_OUT:-$ROOT/output/national_acs_headship_20261112/main}"
export ESTIMATOR_FILE="$ROOT/snapshots/national_acs_f81814dd/code/empirical/acs/kleven_pseudo/estimate_national_first_birth_housing.R"
Rscript --vanilla "$ROOT/code_headship_20261112/headship_estimate.R"
