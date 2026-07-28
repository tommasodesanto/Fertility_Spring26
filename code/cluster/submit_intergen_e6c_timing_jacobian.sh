#!/bin/bash
# Strict five-solve local identification check for the two E6c timing
# coordinates at the certified E6a+E6b winner.
#SBATCH --job-name=ihfe6cjac
#SBATCH --output=logs/slurm_ihfe6cjac_%j.out
#SBATCH --error=logs/slurm_ihfe6cjac_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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

exec "$PYTHON_BIN" \
  "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e6c_timing_jacobian.py" \
  --source "$PROJECT_ROOT/output/model/eqscale_seq_e6ab_recalibration_20260727/report/results.json" \
  --outdir "$PROJECT_ROOT/output/model/eqscale_seq_e6c_timing_jacobian_20260727"
