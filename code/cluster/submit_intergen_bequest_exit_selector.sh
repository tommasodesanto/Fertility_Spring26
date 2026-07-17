#!/bin/bash
#SBATCH --job-name=ihfbqsel
#SBATCH --output=logs/slurm_ihf_bqsel_%j.out
#SBATCH --error=logs/slurm_ihf_bqsel_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"

RESULTS_ROOT="${BEQUEST_EXIT_RESULTS_ROOT:?required}"
WINNER_JSON="${BEQUEST_EXIT_WINNER_JSON:?required}"
exec "$PYTHON_BIN" "$MODEL_DIR/tools/select_intergen_bequest_exit_headline.py" \
  --results-root "$RESULTS_ROOT" \
  --out "$WINNER_JSON" \
  --primary-ltv 0.4 \
  --primary-theta1 0.25
