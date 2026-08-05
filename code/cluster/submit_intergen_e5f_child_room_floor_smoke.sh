#!/bin/bash
# Two-chain E5F Torch smoke. Prepared only; do not submit without authorization.
#SBATCH --job-name=ihfe5fsm
#SBATCH --output=logs/slurm_ihfe5fsm_%A_%a.out
#SBATCH --error=logs/slurm_ihfe5fsm_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:20:00
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
export E3_L4=1 E5=1 E5_MATURATION_REPAIR=1 E5F=1
unset E6A E6B E6C
export E3_TFR_TOP_BIN_WEIGHT="${E3_TFR_TOP_BIN_WEIGHT:-3.602359422009}"
export E2_SEED_RECORD="${E5F_SEED_RECORD:-$PROJECT_ROOT/output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json}"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
RUN_TAG="${E5F_RUN_TAG:-20260805}"
OUTDIR="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_${RUN_TAG}/smoke/chain_${TASK_ID}"
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
  --outdir "$OUTDIR" \
  --seed "$((2026080570 + TASK_ID))" \
  --start-mix "$((TASK_ID - 1))e-1" \
  --minutes 15 \
  --max-evals 10 \
  --smoke \
  --J 17 \
  --Nb 120 \
  --max-iter-eq 40 \
  --tol-eq 2.5e-5
