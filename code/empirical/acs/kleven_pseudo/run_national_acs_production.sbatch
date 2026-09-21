#!/usr/bin/env bash
# Checkpointed national ACS state-by-state matcher, housing bridge, and pool.
# Full state panels remain on disk; the production pool narrows each state to
# the reviewed estimator contract before binding rows.
#SBATCH --job-name=acs-national-production
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_national_production_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_national_production_%j.log
#SBATCH --chdir=/scratch/td2248/projects/kleven_acs_pilot_20260917
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
PARTITION_OUTDIR="${PARTITION_OUTDIR:?completed source-stage output required}"
OUTDIR="${OUTDIR:?fresh national output required}"
case "$PARTITION_OUTDIR" in "$ROOT"/*) ;; *) echo "PARTITION_OUTDIR must be under approved ROOT" >&2; exit 2 ;; esac
case "$OUTDIR" in "$ROOT"/*) ;; *) echo "OUTDIR must be under approved ROOT" >&2; exit 2 ;; esac
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
export PHASE=production RUN_PRODUCTION=YES
Rscript --vanilla "$ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_match_housing.R"
