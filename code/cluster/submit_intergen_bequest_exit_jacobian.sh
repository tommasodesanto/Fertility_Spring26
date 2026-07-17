#!/bin/bash
#SBATCH --job-name=ihfbqjac
#SBATCH --output=logs/slurm_ihf_bqjac_%j.out
#SBATCH --error=logs/slurm_ihf_bqjac_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

WINNER_JSON="${BEQUEST_EXIT_WINNER_JSON:?required}"
OUTDIR="${BEQUEST_EXIT_JACOBIAN_ROOT:?required}"
ARGS=(
  --winner-json "$WINNER_JSON" \
  --outdir "$OUTDIR" \
  --arm "${BEQUEST_EXIT_JACOBIAN_ARM:-A3}" \
  --ltv-terminal 0.4 \
  --theta1 0.25 \
  --unit-step "${BEQUEST_EXIT_JACOBIAN_STEP:-0.005}"
  --max-iter-eq "${BEQUEST_EXIT_MAX_ITER_EQ:-40}"
  --tol-eq "${BEQUEST_EXIT_TOL_EQ:-2.5e-5}"
)
if [ -n "${BEQUEST_EXIT_WINNER_ARM:-}" ]; then
  ARGS+=(--winner-arm "$BEQUEST_EXIT_WINNER_ARM")
fi
if [ "${BEQUEST_EXIT_JACOBIAN_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
if [ -n "${BEQUEST_EXIT_JACOBIAN_PARAMETERS:-}" ]; then
  IFS=',' read -r -a PARAMETERS <<< "$BEQUEST_EXIT_JACOBIAN_PARAMETERS"
  for PARAMETER in "${PARAMETERS[@]}"; do
    ARGS+=(--parameter "$PARAMETER")
  done
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/audit_intergen_bequest_exit_jacobian.py" "${ARGS[@]}"
