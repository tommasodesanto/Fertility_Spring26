#!/bin/bash
# Smoke test for the M5 income-disciplined recalibration: two short --smoke
# chains that exercise the exact production loop structure, one seeded from
# the strict M1 winner and one from the M4 winner, with a nonzero
# tenure_choice_kappa start on the second chain.
#SBATCH --job-name=ihfincdqs
#SBATCH --output=logs/slurm_ihf_incdqs_%A_%a.out
#SBATCH --error=logs/slurm_ihf_incdqs_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:25:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-2

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
RESULTS_ROOT="${INCOME_DISCIPLINED_SMOKE_ROOT:-$RUN_ROOT/smoke}"
M1_SEED_RECORD="${INCOME_DISCIPLINED_M1_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_mortality_recalibration_20260715/report/results.json}"
M4_SEED_RECORD="${INCOME_DISCIPLINED_M4_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_standard_bequest_recalibration_20260716/final_report/results.json}"
ACCEPTANCE_CSV="${INCOME_DISCIPLINED_ACCEPTANCE_CSV:-$PROJECT_ROOT/code/data/mms_center_periphery/output_ownership_audit/acs_ownership_4year_acceptance_bins_6284.csv}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 2 ]; then
  echo "M5 smoke requires SLURM_ARRAY_TASK_ID in 1..2" >&2
  exit 2
fi

CHAIN_IDX=$((TASK_ID - 1))
SEED_ARMS=(M1 M4)
SEED_RECORDS=("$M1_SEED_RECORD" "$M4_SEED_RECORD")
THETA1_STARTS=(0.1 0.4217108770366293)
THETA0_STARTS=(0.0 0.247497728306924)
KAPPA_STARTS=(0.0 0.01)
OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_M5_smoke${CHAIN_IDX}"
mkdir -p "$OUTDIR"

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" \
  --seed-record "${SEED_RECORDS[$CHAIN_IDX]}" \
  --seed-arm "${SEED_ARMS[$CHAIN_IDX]}" \
  --acceptance-csv "$ACCEPTANCE_CSV" \
  --outdir "$OUTDIR" \
  --arm M5 \
  --seed-theta0 "${THETA0_STARTS[$CHAIN_IDX]}" \
  --theta1 "${THETA1_STARTS[$CHAIN_IDX]}" \
  --seed-kappa "${KAPPA_STARTS[$CHAIN_IDX]}" \
  --seed "$((${INCOME_DISCIPLINED_SEED_BASE:-2026071750} + TASK_ID))" \
  --start-mix 0.0 \
  --max-evals 15 \
  --minutes 8 \
  --smoke
