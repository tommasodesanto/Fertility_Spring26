#!/bin/bash
# E6a+E6b: signed twelve-row / ten-parameter combined refit. The production
# seed is the certified E6a winner; smoke uses E5b while E6a is still running.
#SBATCH --job-name=ihfe6ab
#SBATCH --output=logs/slurm_ihfe6ab_%A_%a.out
#SBATCH --error=logs/slurm_ihfe6ab_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
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
export E3_TFR_TOP_BIN_WEIGHT="${E3_TFR_TOP_BIN_WEIGHT:-3.602359422009}"
if [ "${E6AB_SMOKE:-0}" = "1" ]; then
  export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
  RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_smoke_20260727"
else
  export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_e6a_recalibration_20260727/report/results.json"
  RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_recalibration_20260727/production"
fi
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
OUTDIR="$RUN_ROOT/chain_${TASK_ID}"
if [ "$TASK_ID" -le 2 ]; then
  START_MIX=0.0
elif [ "$TASK_ID" -le 5 ]; then
  START_MIX=0.10
else
  START_MIX=0.25
fi
EXTRA_ARGS=()
if [ "${E6AB_SMOKE:-0}" = "1" ]; then
  EXTRA_ARGS=(--smoke --max-evals 13 --minutes 8)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
  --outdir "$OUTDIR" \
  --seed "$((2026072740 + TASK_ID))" \
  --start-mix "$START_MIX" \
  --minutes 225 \
  "${EXTRA_ARGS[@]}"
