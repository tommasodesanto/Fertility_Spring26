#!/bin/bash
# Repaired-E5 funded policies under quota closure, with matched logit
# sensitivities. Six deterministic tasks: three certified arms x two closures.
# Each task runs one baseline plus the two already-approved funded policies.
#SBATCH --job-name=e5quota
#SBATCH --output=logs/slurm_e5quota_%A_%a.out
#SBATCH --error=logs/slurm_e5quota_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-6

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
RUN_ROOT="${E5_POLICY_RUN_ROOT:-$PROJECT_ROOT/output/model/eqscale_seq_e5_policy_quota_closure_20260807/raw_runs}"
case "$TASK_ID" in
  1|4)
    ARM="floor"
    SOURCE="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
    ;;
  2|5)
    ARM="tilt"
    SOURCE="$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_psinneg_extended_20260806/report/results.json"
    ;;
  3|6)
    ARM="chain6"
    SOURCE="$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json"
    ;;
  *)
    echo "unsupported task id: $TASK_ID" >&2
    exit 2
    ;;
esac

if [ "$TASK_ID" -le 3 ]; then
  CLOSURE_MODE="quota"
else
  CLOSURE_MODE="logit"
fi
OUTDIR="$RUN_ROOT/$ARM/$CLOSURE_MODE"
ARGS=(
  --source "$SOURCE"
  --outdir "$OUTDIR"
  --outside-origin-entrant-share 0.169
  --closure-mode "$CLOSURE_MODE"
)
if [ "$CLOSURE_MODE" = "logit" ]; then
  ARGS+=(--entry-taste-scale 2.0)
fi
if [ -n "${E5_POLICY_CASE:-}" ]; then
  ARGS+=(--case "$E5_POLICY_CASE")
fi
if [ "${E5_POLICY_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_e5_repaired_policy_with_entry.py" "${ARGS[@]}"
