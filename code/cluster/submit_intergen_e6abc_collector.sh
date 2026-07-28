#!/bin/bash
# Collect the strict, exactly repeated E6a+E6b+E6c winner.
#SBATCH --job-name=ihfe6abcc
#SBATCH --output=logs/slurm_ihfe6abcc_%j.out
#SBATCH --error=logs/slurm_ihfe6abcc_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"

exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/collect_e1.py" \
  --results-root "$PROJECT_ROOT/output/model/eqscale_seq_e6abc_recalibration_20260727/production" \
  --outdir "$PROJECT_ROOT/output/model/eqscale_seq_e6abc_recalibration_20260727/report"
