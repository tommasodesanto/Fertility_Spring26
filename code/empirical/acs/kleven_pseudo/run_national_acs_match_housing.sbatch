#!/usr/bin/env bash
# State matcher, exact housing bridge, and reviewed pooled-women estimator.
# Smoke uses STATEFIP=50 and performs vanilla-versus-adapter clean/match gates;
# production pools completed state panels after the smoke passes.
#SBATCH --job-name=acs-match-housing
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_match_housing_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_match_housing_%j.log
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
PARTITION_OUTDIR="${PARTITION_OUTDIR:?set to the completed source-stage output}"
OUTDIR="${OUTDIR:?set to a fresh output directory}"
STATEFIP="${STATEFIP:-50}"
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT"
PHASE="${PHASE:-smoke}" STATEFIP="$STATEFIP" PARTITION_OUTDIR="$PARTITION_OUTDIR" OUTDIR="$OUTDIR" Rscript --vanilla code/empirical/acs/kleven_pseudo/run_national_acs_match_housing.R
