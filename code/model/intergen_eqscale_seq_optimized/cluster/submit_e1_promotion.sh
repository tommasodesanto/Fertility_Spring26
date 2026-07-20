#!/usr/bin/env bash
#SBATCH --job-name=e1promote
#SBATCH --output=logs/slurm_e1promote_%A_%a.out
#SBATCH --error=logs/slurm_e1promote_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-1
set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
PROJECT_ROOT="${PROMOTION_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../../.." && pwd)}"
MODEL_DIR="$PROJECT_ROOT/code/model"
TASK="${PROMOTION_TASK:?set PROMOTION_TASK to demand-map, root-case, strict-repeat, or throughput}"
RUN_ROOT="${PROMOTION_RUN_ROOT:-$PROJECT_ROOT/output/model/eqscale_seq_optimized_promotion_20260719}"
mkdir -p "$SCRIPT_DIR/logs" "$RUN_ROOT"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="${PROMOTION_PYTHON:-$(command -v python3 || command -v python)}"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SMOKE_ARGS=(); [ "${PROMOTION_SMOKE:-0}" = 1 ] && SMOKE_ARGS+=(--smoke)
case "$TASK" in
  demand-map|root-case|strict-repeat|throughput)
    exec "$PYTHON_BIN" -m intergen_eqscale_seq_optimized.promotion_worker_e1 "$TASK" --output "$RUN_ROOT/$TASK.json" "${SMOKE_ARGS[@]}"
    ;;
  collect)
    exec "$PYTHON_BIN" -m intergen_eqscale_seq_optimized.promotion_collect_e1 --run-root "$RUN_ROOT"
    ;;
  *) echo "unknown PROMOTION_TASK=$TASK" >&2; exit 2 ;;
esac
