#!/bin/bash
# Strict collector for the August 9 E5b ID-FE refit. Submit afterok on the
# production array; it rejects mixed or stale target contracts.
#SBATCH --job-name=ihfe5bidc
#SBATCH --output=logs/slurm_ihfe5bidc_%j.out
#SBATCH --error=logs/slurm_ihfe5bidc_%j.err
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

RUN_ROOT="${E5B_IDFE_RUN_ROOT:-$PROJECT_ROOT/output/model/eqscale_seq_e5b_idfe_nestingfixed_recalibration_20260809}"
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/collect_e1.py" \
  --results-root "$RUN_ROOT/production" \
  --outdir "$RUN_ROOT/report" \
  --expected-target-set e5_idfe_review_20260809 \
  --expected-target-fingerprint b7c1c1da7578d19e15415c377f2c68b813aca22c6fe5bebda684c14131ec58bc
