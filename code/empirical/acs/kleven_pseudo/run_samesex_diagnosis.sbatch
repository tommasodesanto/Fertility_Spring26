#!/usr/bin/env bash
# SameSex2 negative-rooms diagnostic: ONE gated 4CPU/128GB/<=60min batch
# (no separate 64GB smoke round trip). Runs local fixture tests first,
# then Stage A (identity+numeric reproduction gate), Stage B (event age
# 0:5 RF/FS), Stage C-smoke (one state's exact roster-mapping loop),
# Stage C-full (all 51 states + strict person-key coverage gate), then
# Stage D/E only if all prior gates passed. Uses the UNCHANGED production
# roster/instrument library from the immutable national snapshot
# (twins_samesex_iv_recovery_0f239725) plus this diagnosis pass's own
# reviewed diagnostic files.
#SBATCH --job-name=acs-samesex-diagnosis
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_samesex_diagnosis_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_samesex_diagnosis_%j.log
set -euo pipefail
ROOT_REPO=/scratch/td2248/projects/kleven_acs_pilot_20260917
OUTDIR="${OUTDIR:?set OUTDIR to a fresh directory}"
DIAG_CODE_ROOT="${DIAG_CODE_ROOT:?set DIAG_CODE_ROOT to an immutable diagnosis code snapshot dir}"
NATIONAL_OUTDIR="${NATIONAL_OUTDIR:?set NATIONAL_OUTDIR to the completed national run output dir}"
if [ -d "$OUTDIR" ]; then
  echo "OUTDIR already exists, refusing to reuse: $OUTDIR" >&2
  exit 1
fi
module load r/4.5.1
export R_LIBS_USER="$ROOT_REPO/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT_REPO"

HEARTBEAT_FILE="${OUTDIR%/}_heartbeat_independent.log"
( trap 'echo "$(date -u +%FT%TZ) heartbeat_stopped job=$SLURM_JOB_ID" >> "$HEARTBEAT_FILE"; exit 0' TERM INT EXIT
  while true; do
    echo "$(date -u +%FT%TZ) heartbeat_alive job=$SLURM_JOB_ID outdir=$OUTDIR (monitoring only, not scientific progress)" >> "$HEARTBEAT_FILE"
    sleep 300
  done
) &
HEARTBEAT_PID=$!
cleanup_heartbeat() { kill "$HEARTBEAT_PID" 2>/dev/null || true; }
trap cleanup_heartbeat EXIT

echo "=== stage 0: local fixture tests (meaningful sex/order/link/unequal-weight/group-contrast assertions) ==="
TIMEFORMAT='stage0_tests: elapsed_sec=%R user_sec=%U system_sec=%S'
time { Rscript --vanilla "$DIAG_CODE_ROOT/test_twins_samesex_iv.R" > /dev/null; } 2>&1
time { Rscript --vanilla "$DIAG_CODE_ROOT/test_samesex_diagnosis.R"; } 2>&1

echo "=== stage A-E: run_samesex_diagnosis.R ==="
TIMEFORMAT='samesex_diagnosis_driver: elapsed_sec=%R user_sec=%U system_sec=%S'
time {
ROOT="$DIAG_CODE_ROOT" \
PARTITION_DIR="$ROOT_REPO/output/national_acs_source_stage_20260921/partitions" \
NATIONAL_OUTDIR="$NATIONAL_OUTDIR" \
OUTDIR="$OUTDIR" \
STATUS_PATH="$OUTDIR/status_heartbeat.json" \
STATEFIP_SMOKE=50 \
Rscript --vanilla "$DIAG_CODE_ROOT/run_samesex_diagnosis.R"
} 2>&1

echo "SAMESEX_DIAGNOSIS_BATCH_COMPLETE"
