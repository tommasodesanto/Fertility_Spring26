#!/bin/bash
# E6a+E6b refit under equal-block mean absolute proportional target gaps.
# The model, targets, free parameters, bounds, and strict repeat contract are unchanged.
#SBATCH --job-name=ihfe6abl1
#SBATCH --output=logs/slurm_ihfe6abl1_%A_%a.out
#SBATCH --error=logs/slurm_ihfe6abl1_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-8

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export E3_L4=1 E5=1 E6A=1 E6B=1
unset E6C
export E3_TFR_TOP_BIN_WEIGHT="${E3_TFR_TOP_BIN_WEIGHT:-3.602359422009}"
export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_plainvanilla_local_20260804/report/results.json"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "${E6AB_L1_SMOKE:-0}" = "1" ]; then
  RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_l1_smoke_20260805"
else
  RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_l1_recalibration_20260805/production"
fi
OUTDIR="$RUN_ROOT/chain_${TASK_ID}"
case "$TASK_ID" in
  1) START_MIX=0.00 ;;
  2) START_MIX=0.03 ;;
  3) START_MIX=0.06 ;;
  4) START_MIX=0.10 ;;
  5) START_MIX=0.12 ;;
  6) START_MIX=0.16 ;;
  7) START_MIX=0.20 ;;
  8) START_MIX=0.25 ;;
  *) echo "unsupported task id: $TASK_ID" >&2; exit 2 ;;
esac
EXTRA_ARGS=()
if [ "${E6AB_L1_SMOKE:-0}" = "1" ]; then
  EXTRA_ARGS=(--smoke --max-evals 13 --minutes 8)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
  --outdir "$OUTDIR" \
  --seed "$((2026080500 + TASK_ID))" \
  --start-mix "$START_MIX" \
  --minutes 225 \
  --max-evals 1200 \
  --weight-scheme target_relative_block_equal_l1 \
  "${EXTRA_ARGS[@]}"
