#!/usr/bin/env bash
# Household-position key from extract27 for the headship designs (2026-11-12).
#SBATCH --job-name=acs-headship-key
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_headship_20261112/logs/key_%j.log
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs" OMP_NUM_THREADS=1
export OUTDIR="$ROOT/output/national_acs_headship_20261112/key_build"
Rscript --vanilla "$ROOT/code_headship_20261112/headship_build_key.R"
