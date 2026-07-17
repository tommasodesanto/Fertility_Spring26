#!/bin/bash
#SBATCH --job-name=ihfbqmini
#SBATCH --output=logs/slurm_ihf_bqmini_%j.out
#SBATCH --error=logs/slurm_ihf_bqmini_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RESULTS_ROOT="${BEQUEST_EXERCISE_RESULTS_ROOT:?required}"
SEED_RECORD="${BEQUEST_EXERCISE_SEED_RECORD:?required}"
MODE="${BEQUEST_EXERCISE_MODE:-production}"
mkdir -p "$RESULTS_ROOT"

if [ "$MODE" = smoke ]; then
  exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_calibration_exercise.py" \
    --seed-record "$SEED_RECORD" --outdir "$RESULTS_ROOT" --smoke --quiet
fi

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_calibration_exercise.py" \
  --seed-record "$SEED_RECORD" --outdir "$RESULTS_ROOT" \
  --theta0-grid "${BEQUEST_EXERCISE_THETA0_GRID:-0.05,0.30}" \
  --theta1-grid "${BEQUEST_EXERCISE_THETA1_GRID:-0.25,1.00}" \
  --schedule-arms "${BEQUEST_EXERCISE_SCHEDULE_ARMS:-none,linear_66_82}" \
  --polish-evals "${BEQUEST_EXERCISE_POLISH_EVALS:-12}" \
  --minutes "${BEQUEST_EXERCISE_MINUTES:-35}" --quiet
