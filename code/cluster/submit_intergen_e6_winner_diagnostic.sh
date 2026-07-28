#!/bin/bash
# Strict certified-winner diagnostic packet for E6a, E6b, E6a+E6b, or E6ABC.
# Required environment: E6_DIAGNOSTIC_ARM, E6_DIAGNOSTIC_SOURCE,
# E6_DIAGNOSTIC_OUTDIR.
#SBATCH --job-name=ihfe6diag
#SBATCH --output=logs/slurm_ihfe6diag_%j.out
#SBATCH --error=logs/slurm_ihfe6diag_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:45:00
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

: "${E6_DIAGNOSTIC_ARM:?required: e6a, e6b, e6ab, or e6abc}"
: "${E6_DIAGNOSTIC_SOURCE:?required: collector results.json path}"
: "${E6_DIAGNOSTIC_OUTDIR:?required: output directory}"

exec "$PYTHON_BIN" \
  "$MODEL_DIR/intergen_eqscale_seq_optimized/build_e6_winner_diagnostics.py" \
  --arm "$E6_DIAGNOSTIC_ARM" \
  --source "$E6_DIAGNOSTIC_SOURCE" \
  --outdir "$E6_DIAGNOSTIC_OUTDIR"
