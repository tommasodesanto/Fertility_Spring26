#!/usr/bin/env bash
# Twin1/SameSex2 IV national production run at 128GB/3h -- same 51
# states/sample/science/checkpoints as run_twins_samesex_iv_national.sbatch
# (that 64GB/2h launcher is UNCHANGED and left in place), raised only after
# job 18247804's measured 64GB OOM and the memory/checkpoint fixes
# (274f9c76, 29ae8c77, 698192c8, d4522c51, a7b0c01d). Runs from an
# IMMUTABLE code snapshot and writes an independent 5-minute heartbeat
# outside OUTDIR, monitoring-only.
#SBATCH --job-name=acs-twins-samesex-national-128g
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_national128g_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_national128g_%j.log
set -euo pipefail
ROOT_REPO=/scratch/td2248/projects/kleven_acs_pilot_20260917
OUTDIR="${OUTDIR:?set OUTDIR to a fresh directory}"
CODE_ROOT="${CODE_ROOT:?set CODE_ROOT to an immutable code snapshot dir}"
if [ -d "$OUTDIR" ]; then
  echo "OUTDIR already exists, refusing to reuse: $OUTDIR" >&2
  exit 1
fi
STATEFIP_LIST="1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56"
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

# /usr/bin/time -v is not present on compute nodes here; use bash's
# reserved `time` keyword timing the env-prefixed command group instead
# (wall/user/sys only -- no MaxRSS from this timer). Actual peak memory
# comes from Slurm sacct MaxRSS on the .batch step after the job ends.
TIMEFORMAT='national_driver: elapsed_sec=%R user_sec=%U system_sec=%S'
time {
ROOT="$CODE_ROOT" \
PARTITION_DIR="$ROOT_REPO/output/national_acs_source_stage_20260921/partitions" \
STATEFIP_LIST="$STATEFIP_LIST" \
OUTDIR="$OUTDIR" \
STATUS_PATH="$OUTDIR/status_heartbeat.json" \
SAMPLE_LABEL=young \
Rscript --vanilla "$CODE_ROOT/run_twins_samesex_iv.R"
} 2>&1
