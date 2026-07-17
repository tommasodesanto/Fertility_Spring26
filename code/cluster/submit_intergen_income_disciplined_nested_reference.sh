#!/bin/bash
# Exact theta0=0, tenure_choice_kappa=0 nested M1 reference priced under the
# M5 income-disciplined target system. Theta1 is active but behaviorally
# irrelevant at theta0=0.
#SBATCH --job-name=ihfincdq0
#SBATCH --output=logs/slurm_ihf_incdq0_%j.out
#SBATCH --error=logs/slurm_ihf_incdq0_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general

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
OUTDIR="${INCOME_DISCIPLINED_NESTED_OUTDIR:-$RUN_ROOT/nested_reference}"
M1_SEED_RECORD="${INCOME_DISCIPLINED_M1_SEED_RECORD:-$PROJECT_ROOT/output/model/intergen_mortality_recalibration_20260715/report/results.json}"
ACCEPTANCE_CSV="${INCOME_DISCIPLINED_ACCEPTANCE_CSV:-$PROJECT_ROOT/code/data/mms_center_periphery/output_ownership_audit/acs_ownership_4year_acceptance_bins_6284.csv}"
mkdir -p "$OUTDIR"

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" \
  --seed-record "$M1_SEED_RECORD" \
  --seed-arm M1 \
  --acceptance-csv "$ACCEPTANCE_CSV" \
  --outdir "$OUTDIR" \
  --arm M5 \
  --seed-theta0 0.0 \
  --theta1 0.25 \
  --seed-kappa 0.0 \
  --seed "${INCOME_DISCIPLINED_SEED_BASE:-2026071750}" \
  --start-mix 0.0 \
  --max-evals 1 \
  --minutes 10 \
  --max-iter-eq 10 \
  --tol-eq 1e-4
