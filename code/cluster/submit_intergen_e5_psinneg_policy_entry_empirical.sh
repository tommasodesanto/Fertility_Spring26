#!/bin/bash
# Funded floor-versus-tilt policy comparison using the empirical 16.9 percent
# outside-origin entrant share. The candidate-specific entry probability is
# computed inside the policy driver from each funded baseline flow identity.
#SBATCH --job-name=ihfe5pe
#SBATCH --output=logs/slurm_ihfe5pe_%A_%a.out
#SBATCH --error=logs/slurm_ihfe5pe_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-2

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
RUN_TAG="${E5_POLICY_RUN_TAG:-empirical_entry_20260806}"
case "$TASK_ID" in
  1)
    SOURCE="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
    OUTDIR="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_policy_${RUN_TAG}"
    ;;
  2)
    SOURCE="$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_psinneg_extended_20260806/report/results.json"
    OUTDIR="$PROJECT_ROOT/output/model/eqscale_seq_e5_maturation_repair_psinneg_policy_${RUN_TAG}"
    ;;
  *)
    echo "unsupported task id: $TASK_ID" >&2
    exit 2
    ;;
esac

ARGS=(
  --source "$SOURCE"
  --outdir "$OUTDIR"
  --outside-origin-entrant-share 0.169
  --entry-taste-scale 2.0
)
if [ "${E5_POLICY_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_e5_repaired_policy_with_entry.py" "${ARGS[@]}"
