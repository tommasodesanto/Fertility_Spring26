#!/bin/bash
# Income-disciplined recalibration (arm M5): 11 clean-frontier coordinates plus
# internally estimated theta0, theta1, and tenure_choice_kappa against the
# income-disciplined target set, under the Sommer-Sullivan income anchor and
# the PSID 18-24 entrant wealth distribution. Eight chains use dispersed
# starts with seed arms alternating between the strict M1 and M4 winners.
#SBATCH --job-name=ihfincdq
#SBATCH --output=logs/slurm_ihf_incdq_%A_%a.out
#SBATCH --error=logs/slurm_ihf_incdq_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-8

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RUN_ROOT="$PROJECT_ROOT/output/model/intergen_income_disciplined_recalibration_20260716"
RESULTS_ROOT="${INCOME_DISCIPLINED_RESULTS_ROOT:-$RUN_ROOT/production}"
M1_SEED_RECORD="${INCOME_DISCIPLINED_M1_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_mortality_recalibration_20260715/report/results.json}"
M4_SEED_RECORD="${INCOME_DISCIPLINED_M4_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_standard_bequest_recalibration_20260716/final_report/results.json}"
ACCEPTANCE_CSV="${INCOME_DISCIPLINED_ACCEPTANCE_CSV:-$PROJECT_ROOT/code/data/mms_center_periphery/output_ownership_audit/acs_ownership_4year_acceptance_bins_6284.csv}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 8 ]; then
  echo "M5 requires SLURM_ARRAY_TASK_ID in 1..8" >&2
  exit 2
fi

CHAIN_IDX=$((TASK_ID - 1))
SEED_ARMS=(M1 M4 M1 M4 M1 M4 M1 M4)
SEED_RECORDS=("$M1_SEED_RECORD" "$M4_SEED_RECORD" "$M1_SEED_RECORD" "$M4_SEED_RECORD" "$M1_SEED_RECORD" "$M4_SEED_RECORD" "$M1_SEED_RECORD" "$M4_SEED_RECORD")
THETA1_STARTS=(0.1 0.25 0.4217108770366293 1.0 4.0 12.0 0.25 2.0)
THETA0_STARTS=(0.0 0.1 0.247497728306924 0.5 1.0 2.0 0.247497728306924 0.05)
KAPPA_STARTS=(0.0 0.0 0.01 0.02 0.05 0.0 0.01 0.0)
SEED_ARM="${SEED_ARMS[$CHAIN_IDX]}"
SEED_RECORD="${SEED_RECORDS[$CHAIN_IDX]}"
THETA0_START="${INCOME_DISCIPLINED_SEED_THETA0:-${THETA0_STARTS[$CHAIN_IDX]}}"
THETA1_START="${INCOME_DISCIPLINED_SEED_THETA1:-${THETA1_STARTS[$CHAIN_IDX]}}"
KAPPA_START="${INCOME_DISCIPLINED_SEED_KAPPA:-${KAPPA_STARTS[$CHAIN_IDX]}}"
OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_M5_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"

ARGS=(
  --seed-record "$SEED_RECORD"
  --seed-arm "$SEED_ARM"
  --acceptance-csv "$ACCEPTANCE_CSV"
  --outdir "$OUTDIR"
  --arm M5
  --seed-theta0 "$THETA0_START"
  --theta1 "$THETA1_START"
  --seed-kappa "$KAPPA_START"
  --seed "$((${INCOME_DISCIPLINED_SEED_BASE:-2026071750} + TASK_ID))"
  --start-mix "${INCOME_DISCIPLINED_START_MIX:-0.00}"
  --max-evals "${INCOME_DISCIPLINED_MAX_EVALS:-1000}"
  --minutes "${INCOME_DISCIPLINED_MINUTES:-225}"
  --max-iter-eq "${INCOME_DISCIPLINED_MAX_ITER_EQ:-10}"
  --tol-eq "${INCOME_DISCIPLINED_TOL_EQ:-1e-4}"
)
if [ "${INCOME_DISCIPLINED_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" "${ARGS[@]}"
