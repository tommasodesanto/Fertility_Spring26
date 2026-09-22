#!/usr/bin/env bash
# Estimator-only national ACS continuation.  This launcher must be reviewed
# after the production matcher is terminal; it never invokes matching or raw
# ACS/CPS preparation.
#SBATCH --job-name=acs-national-estimator-cont
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_national_estimator_cont_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_national_estimator_cont_%j.log
#SBATCH --chdir=/scratch/td2248/projects/kleven_acs_pilot_20260917
set -euo pipefail

ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
SNAPSHOT_ROOT="${SNAPSHOT_ROOT:?immutable reviewed ACS snapshot required}"
CONTINUATION_CODE_ROOT="${CONTINUATION_CODE_ROOT:-$SNAPSHOT_ROOT}"
PRODUCTION_DIR="${PRODUCTION_DIR:-$ROOT/output/national_acs_production_20260921a}"
OUTDIR="${OUTDIR:-$ROOT/output/national_acs_estimator_continuation_20260921a}"
PANEL_FILE="${PANEL_FILE:-$PRODUCTION_DIR/national_cps_acs_pseudo-panel_housing.rds}"
REUSE_DIR="${REUSE_DIR:-$PRODUCTION_DIR/national_first_birth_housing}"

case "$SNAPSHOT_ROOT" in "$ROOT"/*) ;; *) echo "SNAPSHOT_ROOT must be under approved ROOT" >&2; exit 2 ;; esac
case "$CONTINUATION_CODE_ROOT" in "$ROOT"/*) ;; *) echo "CONTINUATION_CODE_ROOT must be under approved ROOT" >&2; exit 2 ;; esac
case "$PRODUCTION_DIR" in "$ROOT"/*) ;; *) echo "PRODUCTION_DIR must be under approved ROOT" >&2; exit 2 ;; esac
case "$OUTDIR" in "$ROOT"/*) ;; *) echo "OUTDIR must be under approved ROOT" >&2; exit 2 ;; esac
test -r "$CONTINUATION_CODE_ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_estimator_continuation.R"
test -r "$CONTINUATION_CODE_ROOT/code/empirical/acs/kleven_pseudo/estimate_national_first_birth_housing.R"
test -r "$SNAPSHOT_ROOT/code/empirical/acs/kleven_pseudo/source_contract.json"
test -r "$PANEL_FILE"
test -r "$PRODUCTION_DIR/production_start_receipt.json"
test -r "$PRODUCTION_DIR/national_pool_manifest.json"
test -r "$PRODUCTION_DIR/stage_receipt.json"
test -d "$REUSE_DIR"

module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
export PROJECT_ROOT="$ROOT" SNAPSHOT_ROOT="$SNAPSHOT_ROOT" CONTINUATION_CODE_ROOT="$CONTINUATION_CODE_ROOT"
export PRODUCTION_DIR="$PRODUCTION_DIR" OUTDIR="$OUTDIR" PANEL_FILE="$PANEL_FILE" REUSE_DIR="$REUSE_DIR"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

Rscript --vanilla "$CONTINUATION_CODE_ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_estimator_continuation.R"
