#!/usr/bin/env bash
# Array over 51 states: household-position filter and 1/k' reweighting (2026-11-12).
#SBATCH --job-name=acs-headship-filter
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --array=0-50%17
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_headship_20261112/logs/filter_%A_%a.log
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
FIPS=(1 2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56)
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs" OMP_NUM_THREADS=1
export STATEFIP=${FIPS[$SLURM_ARRAY_TASK_ID]}
export KEY_DIR="$ROOT/output/national_acs_headship_20261112/key_local_mac/key"
export OUTDIR="${HEADSHIP_OUT:-$ROOT/output/national_acs_headship_20261112/main}"
export HH_GQ="${HH_GQ:?set household GQ codes}"
Rscript --vanilla "$ROOT/code_headship_20261112/headship_filter_reweight.R"
