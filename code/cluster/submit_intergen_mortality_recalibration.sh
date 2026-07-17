#!/bin/bash
#SBATCH --job-name=ihfmort
#SBATCH --output=logs/slurm_ihf_mort_%A_%a.out
#SBATCH --error=logs/slurm_ihf_mort_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:40:00
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

RESULTS_ROOT="${MORTALITY_RESULTS_ROOT:?required}"
SEED_RECORD="${MORTALITY_SEED_RECORD:?required}"
ACCEPTANCE_CSV="${MORTALITY_ACCEPTANCE_CSV:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"

# Four fully re-estimated chains per arm.  M0 is the inert no-mortality control;
# M1 adds only the externally pinned SSA post-retirement survival schedule.
CHAIN_IDX=$(((TASK_ID - 1) % 4))
ARM=M0
[ "$TASK_ID" -le 4 ] || ARM=M1
START_MIXES=(0.00 0.04 0.08 0.12)
START_MIX="${START_MIXES[$CHAIN_IDX]}"

OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_${ARM}_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"
ARGS=(
  --seed-record "$SEED_RECORD"
  --acceptance-csv "$ACCEPTANCE_CSV"
  --outdir "$OUTDIR"
  --arm "$ARM"
  --seed "$((${MORTALITY_SEED_BASE:-2026071500} + TASK_ID))"
  --start-mix "$START_MIX"
  --max-evals "${MORTALITY_MAX_EVALS:-1000}"
  --minutes "${MORTALITY_MINUTES:-90}"
  --max-iter-eq "${MORTALITY_MAX_ITER_EQ:-10}"
  --tol-eq "${MORTALITY_TOL_EQ:-1e-4}"
)
if [ "${MORTALITY_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" "${ARGS[@]}"
