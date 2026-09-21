#!/usr/bin/env bash
# Checkpointed national ACS state-by-state matcher, housing bridge, and pool.
# Full state panels remain on disk; the production pool narrows each state to
# the reviewed estimator contract before binding rows.
#SBATCH --job-name=acs-national-production
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_national_production_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_national_production_%j.log
#SBATCH --chdir=/scratch/td2248/projects/kleven_acs_pilot_20260917
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
SNAPSHOT_ROOT="${SNAPSHOT_ROOT:?immutable ACS snapshot required}"
PARTITION_OUTDIR="${PARTITION_OUTDIR:?completed source-stage output required}"
OUTDIR="${OUTDIR:?fresh national output required}"
ALL_STATES="1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56"
case "$PARTITION_OUTDIR" in "$ROOT"/*) ;; *) echo "PARTITION_OUTDIR must be under approved ROOT" >&2; exit 2 ;; esac
case "$OUTDIR" in "$ROOT"/*) ;; *) echo "OUTDIR must be under approved ROOT" >&2; exit 2 ;; esac
case "$SNAPSHOT_ROOT" in "$ROOT"/*) ;; *) echo "SNAPSHOT_ROOT must be under approved ROOT" >&2; exit 2 ;; esac
test -r "$SNAPSHOT_ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_match_housing.R"
test -r "$SNAPSHOT_ROOT/code/empirical/acs/kleven_pseudo/estimate_national_first_birth_housing.R"
test -r "$SNAPSHOT_ROOT/code/empirical/acs/kleven_pseudo/source_contract.json"
test -r "$SNAPSHOT_ROOT/vendor/clean_acs.R"
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
export PROJECT_ROOT="$SNAPSHOT_ROOT" VENDOR_DIR="$SNAPSHOT_ROOT/vendor"
export SOURCE_CONTRACT="$SNAPSHOT_ROOT/code/empirical/acs/kleven_pseudo/source_contract.json"
export HOUSING_RAW="$ROOT/inputs/ACS/local_extract27_20260919/extract27.dta"
export PHASE=production RUN_PRODUCTION=YES STATEFIP="$ALL_STATES"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
HEARTBEAT_FILE="$ROOT/output/acs_national_production_${SLURM_JOB_ID:-manual}.heartbeat"
Rscript --vanilla "$SNAPSHOT_ROOT/code/empirical/acs/kleven_pseudo/run_national_acs_match_housing.R" &
PID=$!
(
  while kill -0 "$PID" 2>/dev/null; do
    printf '%s production_heartbeat pid=%s outdir=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$PID" "$OUTDIR" >> "$HEARTBEAT_FILE"
    sleep 300
  done
) &
HEARTBEAT=$!
set +e
wait "$PID"
STATUS=$?
set -e
kill "$HEARTBEAT" 2>/dev/null || true
wait "$HEARTBEAT" 2>/dev/null || true
exit "$STATUS"
