#!/usr/bin/env bash
# Estimator-only recovery from a saved state panel; no raw read or matching.
# The panel and output paths are supplied as absolute paths under ROOT.
#SBATCH --job-name=acs-est-vt-recovery
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_est_vt_recovery_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_est_vt_recovery_%j.log
#SBATCH --chdir=/scratch/td2248/projects/kleven_acs_pilot_20260917
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
PANEL_FILE="${PANEL_FILE:?absolute saved panel required}"
OUTDIR="${OUTDIR:?absolute fresh output required}"
case "$PANEL_FILE" in "$ROOT"/*) ;; *) echo "PANEL_FILE must be under approved ROOT" >&2; exit 2 ;; esac
case "$OUTDIR" in "$ROOT"/*) ;; *) echo "OUTDIR must be under approved ROOT" >&2; exit 2 ;; esac
test -r "$PANEL_FILE"
test ! -e "$OUTDIR"
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
test -r "$ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_estimator_recovery.R"
Rscript --vanilla "$ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_estimator_recovery.R"
