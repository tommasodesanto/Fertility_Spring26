#!/bin/bash
# Refit the maintained pre-E6 E5b architecture under the August 9 ID-FE
# first-birth rooms target contract.  Production: eight independent local
# chains, at most 1,000 solves or 225 minutes each.  Every solve checkpoints.
# Submit the two-chain smoke and one-chain exact preflight before production.
#SBATCH --job-name=ihfe5bid
#SBATCH --output=logs/slurm_ihfe5bid_%A_%a.out
#SBATCH --error=logs/slurm_ihfe5bid_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-8

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

# Exact maintained architecture: historical E5b, with all later experimental
# model switches forced off.  The target-set name and fingerprint are hard
# gates, so target or weight drift fails before the first model solve.
export E3_L4=1
export E5=1
unset E5_MATURATION_REPAIR E5F E6A E6B E6C E5_PSI_MIN E5_PSI_BOUND || true
export E3_TFR_TOP_BIN_WEIGHT=3.602359422009
export E5_EXPECTED_TARGET_SET=e5_idfe_review_20260809
export E5_EXPECTED_TARGET_FINGERPRINT=b7c1c1da7578d19e15415c377f2c68b813aca22c6fe5bebda684c14131ec58bc
export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
MODE="${E5B_IDFE_MODE:-production}"
RUN_ROOT="${E5B_IDFE_RUN_ROOT:-$PROJECT_ROOT/output/model/eqscale_seq_e5b_idfe_nestingfixed_recalibration_20260809}"

case "$TASK_ID" in
  1) START_MIX=0.000 ;;
  2|3) START_MIX=0.010 ;;
  4|5) START_MIX=0.025 ;;
  *) START_MIX=0.050 ;;
esac

if [ "$MODE" = "smoke" ]; then
  OUTDIR="$RUN_ROOT/smoke/chain_${TASK_ID}"
  EXTRA_ARGS=(--smoke --minutes 8 --max-evals 13)
elif [ "$MODE" = "preflight" ]; then
  if [ "$TASK_ID" != "1" ]; then
    echo "E5b ID-FE preflight must use array task 1 only" >&2
    exit 2
  fi
  OUTDIR="$RUN_ROOT/preflight/chain_1"
  START_MIX=0.000
  EXTRA_ARGS=(
    --minutes 5 --max-evals 1
    --seed-reproduction-record "$E2_SEED_RECORD"
    --seed-reproduction-atol 2e-4
  )
elif [ "$MODE" = "production" ]; then
  OUTDIR="$RUN_ROOT/production/chain_${TASK_ID}"
  EXTRA_ARGS=(--minutes 225 --max-evals 1000)
else
  echo "Unknown E5B_IDFE_MODE=$MODE" >&2
  exit 2
fi

exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
  --outdir "$OUTDIR" \
  --seed "$((2026080900 + TASK_ID))" \
  --start-mix "$START_MIX" \
  "${EXTRA_ARGS[@]}"
