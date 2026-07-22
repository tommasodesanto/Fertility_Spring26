#!/bin/bash
# E3b: overnight continuation of the E3 literal-parity recalibration.
# Eight fresh NM chains seeded from the collected E3 winner (report written
# by the dependent collector); same L4 conventions and 15-moment system.
#SBATCH --job-name=ihfe3b
#SBATCH --output=logs/slurm_ihfe3b_%A_%a.out
#SBATCH --error=logs/slurm_ihfe3b_%A_%a.err
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
export E3_L4=1
export E3_TFR_TOP_BIN_WEIGHT="${E3_TFR_TOP_BIN_WEIGHT:-3.4}"
export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_l4_recalibration_20260722/report/results.json"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_l4_recalibration_20260722/continuation"
OUTDIR="$RUN_ROOT/chain_${TASK_ID}"
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
  --outdir "$OUTDIR" \
  --seed "$((2026072300 + TASK_ID))" \
  --minutes 225
