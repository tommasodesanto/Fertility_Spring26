#!/bin/bash
# Re-estimate the unchanged E5 twelve-target, ten-parameter contract after the
# author-directed independent child-maturation repair.
#SBATCH --job-name=ihfe5mr
#SBATCH --output=logs/slurm_ihfe5mr_%A_%a.out
#SBATCH --error=logs/slurm_ihfe5mr_%A_%a.err
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
export E3_L4=1 E5=1 E5_MATURATION_REPAIR=1
export E3_TFR_TOP_BIN_WEIGHT="${E3_TFR_TOP_BIN_WEIGHT:-3.602359422009}"
export E2_SEED_RECORD="${E5_REPAIR_SEED_RECORD:-$PROJECT_ROOT/output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json}"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
RUN_TAG="${E5_REPAIR_RUN_TAG:-20260805}"
RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_${RUN_TAG}/production"
OUTDIR="$RUN_ROOT/chain_${TASK_ID}"
ARGS=(
  --outdir "$OUTDIR"
  --seed "$((2026080500 + TASK_ID))"
  --minutes "${E5_REPAIR_MINUTES:-225}"
  --max-evals "${E5_REPAIR_MAX_EVALS:-1000}"
)
if [ "${E5_REPAIR_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" "${ARGS[@]}"
