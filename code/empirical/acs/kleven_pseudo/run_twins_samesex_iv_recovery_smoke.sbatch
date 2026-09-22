#!/usr/bin/env bash
# Recovery smoke after the 18247804 OOM + memory/checkpoint fixes
# (274f9c76, 29ae8c77, 698192c8, d4522c51, a7b0c01d): runs the real VT
# driver end to end, then a compute-only synthetic memory probe that
# replicates VT's saved slim analytic frames up to realistic national
# Twin1 scale (857,607 rows) to stress-test the SAME production
# fit_instrument_outcome()/ar_confidence_set() code path. NOT a national
# estimate -- a measured go/no-go gate before the 128GB national attempt.
#SBATCH --job-name=acs-twins-samesex-recovery-smoke
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_recovery_smoke_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_recovery_smoke_%j.log
set -euo pipefail
ROOT_REPO=/scratch/td2248/projects/kleven_acs_pilot_20260917
OUTDIR="${OUTDIR:?set OUTDIR to a fresh directory}"
CODE_ROOT="${CODE_ROOT:?set CODE_ROOT to an immutable code snapshot dir}"
module load r/4.5.1
export R_LIBS_USER="$ROOT_REPO/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT_REPO"

VT_OUTDIR="$OUTDIR/vt_smoke"
PROBE_OUTDIR="$OUTDIR/scaled_probe"
mkdir -p "$VT_OUTDIR" "$PROBE_OUTDIR"

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

echo "=== stage 1: real VT driver (STATEFIP_LIST=50) ==="
ROOT="$CODE_ROOT" \
PARTITION_DIR="$ROOT_REPO/output/national_acs_source_stage_20260921/partitions" \
STATEFIP_LIST=50 \
OUTDIR="$VT_OUTDIR" \
STATUS_PATH="$VT_OUTDIR/status_heartbeat.json" \
SAMPLE_LABEL=young \
/usr/bin/time -v Rscript --vanilla "$CODE_ROOT/run_twins_samesex_iv.R" 2>&1

echo "=== gate: VT must show 18 full_fit + 18 primary_receipts + analytic identity ==="
python3 - "$VT_OUTDIR" << 'PYEOF'
import json, sys, os
d = sys.argv[1]
tally = json.load(open(os.path.join(d, "fit_status_tally.json")))
full = sum(r["N"] for r in tally if r["status"] == "full_fit")
assert full == 18, f"expected 18 full_fit, got {full}: {tally}"
n_primary = len([f for f in os.listdir(os.path.join(d, "primary_receipts")) if f.endswith(".json")])
assert n_primary == 18, f"expected 18 primary_receipts, got {n_primary}"
assert os.path.exists(os.path.join(d, "analytic_frames_identity.json")), "missing analytic_frames_identity.json"
assert os.path.exists(os.path.join(d, "analytic_frames.rds")), "missing analytic_frames.rds"
print("VT_GATE_PASS: 18/18 full_fit, 18 primary_receipts, analytic frames+identity present")
PYEOF

echo "=== stage 2: compute-only scaled memory probe (target Twin1 N=857607) ==="
ROOT="$CODE_ROOT" \
VT_OUTDIR="$VT_OUTDIR" \
PROBE_OUTDIR="$PROBE_OUTDIR" \
TARGET_TWIN1_N=857607 \
/usr/bin/time -v Rscript --vanilla "$CODE_ROOT/probe_scaled_memory.R" 2>&1

echo "RECOVERY_SMOKE_COMPLETE"
