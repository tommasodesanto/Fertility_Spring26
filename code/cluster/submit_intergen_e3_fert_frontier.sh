#!/bin/bash
# Fertility-frontier reachability probe: 147 fixed-coordinate GE cells over
# (psi_child x kappa_fert x gamma_e) at the certified E3 winner, L4 on.
#SBATCH --job-name=ihfe3ff
#SBATCH --output=logs/slurm_ihfe3ff_%A_%a.out
#SBATCH --error=logs/slurm_ihfe3ff_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-147

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_fertility_frontier_scan.py" \
  --cell "${SLURM_ARRAY_TASK_ID:?}" \
  --outdir "$PROJECT_ROOT/output/model/eqscale_fertility_frontier_20260722"
