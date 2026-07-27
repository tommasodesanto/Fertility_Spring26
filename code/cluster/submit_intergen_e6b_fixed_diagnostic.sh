#!/bin/bash
# Fixed-E5b strict comparison of E6b alone and E6a+E6b.
#SBATCH --job-name=ihfe6bd
#SBATCH --output=logs/slurm_ihfe6bd_%j.out
#SBATCH --error=logs/slurm_ihfe6bd_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
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

EXTRA_ARGS=()
OUTDIR="$PROJECT_ROOT/output/model/eqscale_seq_e6b_income_level_diagnostic_20260727"
if [ "${E6B_DIAGNOSTIC_SMOKE:-0}" = "1" ]; then
  EXTRA_ARGS=(--smoke)
  OUTDIR="$PROJECT_ROOT/output/model/eqscale_seq_e6b_income_level_diagnostic_smoke_20260727"
fi

exec "$PYTHON_BIN" \
  "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e6b_income_level_diagnostic.py" \
  --outdir "$OUTDIR" \
  "${EXTRA_ARGS[@]}"
