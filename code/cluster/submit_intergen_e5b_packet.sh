#!/bin/bash
# E5b policy packet: standard diagnostics + policy arms at the certified E5b
# winner (verification gate at 1e-6 against the certified target-fit table).
#SBATCH --job-name=ihfe5bp
#SBATCH --output=logs/slurm_ihfe5bp_%j.out
#SBATCH --error=logs/slurm_ihfe5bp_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/build_e2_packet.py" --arm e5
