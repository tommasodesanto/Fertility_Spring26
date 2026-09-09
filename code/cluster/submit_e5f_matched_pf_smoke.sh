#!/usr/bin/env bash
#SBATCH --job-name=e5f_matched_pf_smoke
#SBATCH --output=logs/e5f_matched_pf_smoke_%A_%a.out
#SBATCH --error=logs/e5f_matched_pf_smoke_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-1%2

# Invoke from code/cluster after creating its logs directory. The contract must
# contain cluster-absolute input paths and pins from this exact source snapshot.
set -euo pipefail
: "${SLURM_JOB_ID:?Submit this file with sbatch}"
: "${SLURM_ARRAY_TASK_ID:?Expected two-arm array}"
: "${E5F_MATCHED_PF_CONTRACT:?Absolute contract path required}"
: "${E5F_MATCHED_PF_CONTRACT_SHA256:?Contract SHA-256 required}"
: "${E5F_MATCHED_PF_OUTPUT:?Absolute new output root required}"
: "${E5F_MATCHED_PF_ROOT:?Absolute source snapshot root required}"
case "$SLURM_ARRAY_TASK_ID" in
    0) ARM=sequential ;;
    1) ARM=nested ;;
    *) exit 2 ;;
esac
for value in "$E5F_MATCHED_PF_CONTRACT" "$E5F_MATCHED_PF_OUTPUT" "$E5F_MATCHED_PF_ROOT"; do
    [[ "$value" = /* ]] || { echo 'All paths must be absolute' >&2; exit 2; }
done
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="$E5F_MATCHED_PF_ROOT/output/cache/matplotlib/$ARM"
export NUMBA_CACHE_DIR="$E5F_MATCHED_PF_ROOT/output/cache/numba/$ARM"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
cd "$E5F_MATCHED_PF_ROOT"
exec "$PYTHON_BIN" code/model/tools/run_e5f_matched_pf_smoke.py \
    --contract "$E5F_MATCHED_PF_CONTRACT" \
    --contract-sha256 "$E5F_MATCHED_PF_CONTRACT_SHA256" \
    --arm "$ARM" --output "$E5F_MATCHED_PF_OUTPUT/$ARM"
