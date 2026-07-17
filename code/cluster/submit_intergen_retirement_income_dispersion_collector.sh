#!/bin/bash
#SBATCH --job-name=ihfretincc
#SBATCH --output=logs/slurm_ihf_retincc_%j.out
#SBATCH --error=logs/slurm_ihf_retincc_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"

RESULTS_ROOT="${RETIREMENT_DISP_RESULTS_ROOT:?required}"
REPORT_DIR="${RETIREMENT_DISP_REPORT_DIR:?required}"
exec "$PYTHON_BIN" "$MODEL_DIR/tools/collect_intergen_retirement_income_dispersion.py" \
  --results-dir "$RESULTS_ROOT" \
  --outdir "$REPORT_DIR"
