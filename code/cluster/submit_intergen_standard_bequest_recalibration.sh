#!/bin/bash
# Standard bequest recalibration (arm M4): 11 clean-frontier coordinates plus
# internally estimated theta0 and theta1 against the median-composition target
# set. Six chains use dispersed starts; theta1=0.25 is only one start.
#SBATCH --job-name=ihfstdbq
#SBATCH --output=logs/slurm_ihf_stdbq_%A_%a.out
#SBATCH --error=logs/slurm_ihf_stdbq_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:25:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-6

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RUN_ROOT="$PROJECT_ROOT/output/model/intergen_standard_bequest_recalibration_20260716"
RESULTS_ROOT="${STANDARD_BEQUEST_RESULTS_ROOT:-$RUN_ROOT/production}"
M1_SEED_RECORD="${STANDARD_BEQUEST_M1_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_mortality_recalibration_20260715/report/results.json}"
M3_SEED_RECORD="${STANDARD_BEQUEST_M3_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_internal_bequest_recalibration_20260715/report/results.json}"
ACCEPTANCE_CSV="${STANDARD_BEQUEST_ACCEPTANCE_CSV:-$PROJECT_ROOT/code/data/mms_center_periphery/output_ownership_audit/acs_ownership_4year_acceptance_bins_6284.csv}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 6 ]; then
  echo "M4 requires SLURM_ARRAY_TASK_ID in 1..6" >&2
  exit 2
fi

CHAIN_IDX=$((TASK_ID - 1))
START_MIXES=(0.00 0.02 0.04 0.06 0.00 0.04)
SEED_ARMS=(M1 M1 M1 M1 M3 M3)
SEED_RECORDS=("$M1_SEED_RECORD" "$M1_SEED_RECORD" "$M1_SEED_RECORD" "$M1_SEED_RECORD" "$M3_SEED_RECORD" "$M3_SEED_RECORD")
THETA0_STARTS=(0.0 0.05 0.20 0.75 0.31024680236397595 2.0)
THETA1_STARTS=(0.25 0.25 1.0 4.0 0.5361411832232645 12.0)
START_MIX="${STANDARD_BEQUEST_START_MIX:-${START_MIXES[$CHAIN_IDX]}}"
SEED_ARM="${SEED_ARMS[$CHAIN_IDX]}"
SEED_RECORD="${SEED_RECORDS[$CHAIN_IDX]}"
THETA0_START="${STANDARD_BEQUEST_SEED_THETA0:-${THETA0_STARTS[$CHAIN_IDX]}}"
THETA1_START="${STANDARD_BEQUEST_SEED_THETA1:-${THETA1_STARTS[$CHAIN_IDX]}}"
OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_M4_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"

ARGS=(
  --seed-record "$SEED_RECORD"
  --seed-arm "$SEED_ARM"
  --acceptance-csv "$ACCEPTANCE_CSV"
  --outdir "$OUTDIR"
  --arm M4
  --seed-theta0 "$THETA0_START"
  --theta1 "$THETA1_START"
  --seed "$((${STANDARD_BEQUEST_SEED_BASE:-2026071640} + TASK_ID))"
  --start-mix "$START_MIX"
  --max-evals "${STANDARD_BEQUEST_MAX_EVALS:-1000}"
  --minutes "${STANDARD_BEQUEST_MINUTES:-75}"
  --max-iter-eq "${STANDARD_BEQUEST_MAX_ITER_EQ:-10}"
  --tol-eq "${STANDARD_BEQUEST_TOL_EQ:-1e-4}"
)
if [ "${STANDARD_BEQUEST_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" "${ARGS[@]}"
