#!/bin/bash
# E5 fixed-kappa_fert probe, with psi_child widened to [-6, 6].
#SBATCH --job-name=ihfe5p
#SBATCH --output=logs/slurm_ihfe5p_%A_%a.out
#SBATCH --error=logs/slurm_ihfe5p_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-10

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
KE_VALUES=(0.5 1.0 2.0 4.0 8.0 12.0 18.0 24.0 30.0 36.339009750351664)
if (( TASK_ID < 1 || TASK_ID > ${#KE_VALUES[@]} )); then
  echo "Expected SLURM_ARRAY_TASK_ID in 1..${#KE_VALUES[@]}" >&2
  exit 2
fi
export E3_L4=1 E5=1 E3_TFR_TOP_BIN_WEIGHT=3.602359422009
export E5_PROBE_FIX_KE="${KE_VALUES[$((TASK_ID - 1))]}" E5_PSI_BOUND=6.0
export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_e5_recalibration_20260724/report/results.json"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e5_probe_20260725/production"
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
  --outdir "$RUN_ROOT/chain_${TASK_ID}" \
  --seed "$((2026072510 + TASK_ID))" \
  --minutes 200
