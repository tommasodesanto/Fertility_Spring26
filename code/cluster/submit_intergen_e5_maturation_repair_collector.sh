#!/bin/bash
# Strict collector for the independent child-maturation E5 repair.
#SBATCH --job-name=ihfe5mrc
#SBATCH --output=logs/slurm_ihfe5mrc_%A.out
#SBATCH --error=logs/slurm_ihfe5mrc_%A.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
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

RUN_TAG="${E5_REPAIR_RUN_TAG:-20260805}"
RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_${RUN_TAG}/production"
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/collect_e1.py" \
  --results-root "$RUN_ROOT" \
  --outdir "$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_${RUN_TAG}/report"
