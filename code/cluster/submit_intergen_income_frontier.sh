#!/bin/bash
# Fixed-M5 income-risk, transfer-floor, and subsistence feasibility frontier.
# One cpu_short task evaluates the 28 cells sequentially; FRONTIER_CELLS can
# restrict the run to a comma-separated smoke-test subset.
#SBATCH --job-name=ihffront
#SBATCH --output=logs/slurm_ihf_front_%j.out
#SBATCH --error=logs/slurm_ihf_front_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
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

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_income_feasibility_frontier.py" ${FRONTIER_CELLS:+--cells "$FRONTIER_CELLS"} --outdir "$PROJECT_ROOT/output/model/income_feasibility_frontier_20260718"
