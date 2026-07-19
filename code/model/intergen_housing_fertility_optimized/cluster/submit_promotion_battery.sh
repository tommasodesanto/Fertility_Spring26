#!/bin/bash
# Isolated optimized-package promotion battery. Select the task through
# PROMOTION_TASK and override --array at sbatch submission.
#SBATCH --job-name=ihfoptgate
#SBATCH --output=logs/slurm_ihfoptgate_%A_%a.out
#SBATCH --error=logs/slurm_ihfoptgate_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-1

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
PROJECT_ROOT="${PROMOTION_PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../../../.." && pwd)}"
MODEL_DIR="$PROJECT_ROOT/code/model"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
TASK="${PROMOTION_TASK:?set PROMOTION_TASK}"
RUN_ROOT="${PROMOTION_RUN_ROOT:-$PROJECT_ROOT/output/model/intergen_optimized_promotion_20260719}"

mkdir -p "$SCRIPT_DIR/logs" "$RUN_ROOT"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SMOKE_ARGS=()
if [ "${PROMOTION_SMOKE:-0}" = "1" ]; then
  SMOKE_ARGS+=(--smoke)
fi

case "$TASK" in
  bound-parity)
    OUTFILE="$RUN_ROOT/bounds/task_${TASK_ID}.json"
    ARGS=(bound-parity --index "$TASK_ID" --output "$OUTFILE")
    ;;
  demand-point)
    OUTFILE="$RUN_ROOT/demand/task_${TASK_ID}.json"
    ARGS=(demand-point --index "$TASK_ID" --output "$OUTFILE")
    ;;
  root-case)
    OUTFILE="$RUN_ROOT/roots/task_${TASK_ID}.json"
    ARGS=(root-case --index "$TASK_ID" --output "$OUTFILE")
    ;;
  strict-repeat)
    OUTFILE="$RUN_ROOT/repeats/task_${TASK_ID}.json"
    ARGS=(strict-repeat --index "$TASK_ID" --output "$OUTFILE")
    ;;
  calibration)
    if [ "$TASK_ID" -eq 1 ]; then
      PACKAGE=production
    elif [ "$TASK_ID" -eq 2 ]; then
      PACKAGE=optimized
    else
      echo "calibration task requires array index 1 or 2" >&2
      exit 2
    fi
    OUTFILE="$RUN_ROOT/calibration/${PACKAGE}.json"
    ARGS=(calibration --package "$PACKAGE" --output "$OUTFILE")
    ;;
  diagnostics)
    exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.diagnostic_parity \
      --outdir "$RUN_ROOT/diagnostics" "${SMOKE_ARGS[@]}"
    ;;
  *)
    echo "unknown PROMOTION_TASK=$TASK" >&2
    exit 2
    ;;
esac

exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.promotion_worker \
  "${ARGS[@]}" "${SMOKE_ARGS[@]}"
