#!/bin/bash
# Strict collector for the M5 income-disciplined recalibration. Fails unless a
# strict, exactly repeated winner exists and it beats the strict M1 theta0=0,
# tenure_choice_kappa=0 nested seed under the M5 objective. The
# identification_reported acceptance row is marked pending when the Jacobian
# summary has not been produced yet.
#SBATCH --job-name=ihfincdqc
#SBATCH --output=logs/slurm_ihf_incdqc_%j.out
#SBATCH --error=logs/slurm_ihf_incdqc_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"

RUN_ROOT="$PROJECT_ROOT/output/model/intergen_income_disciplined_recalibration_20260716"
RESULTS_ROOT="${INCOME_DISCIPLINED_RESULTS_ROOT:-$RUN_ROOT/production}"
REPORT_DIR="${INCOME_DISCIPLINED_REPORT_DIR:-$RUN_ROOT/report}"
NESTED_REFERENCE="${INCOME_DISCIPLINED_NESTED_REFERENCE:-$RUN_ROOT/nested_reference/summary.json}"
JACOBIAN_SUMMARY="${INCOME_DISCIPLINED_JACOBIAN_SUMMARY:-$RUN_ROOT/identification/summary.json}"
exec "$PYTHON_BIN" "$MODEL_DIR/tools/collect_intergen_internal_bequest_recalibration.py" \
  --arm M5 \
  --results-dir "$RESULTS_ROOT" \
  --nested-reference "$NESTED_REFERENCE" \
  --jacobian-summary "$JACOBIAN_SUMMARY" \
  --outdir "$REPORT_DIR"
