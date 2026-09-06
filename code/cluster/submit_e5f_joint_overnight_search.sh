#!/bin/bash
#SBATCH --job-name=e5f_joint_night
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=08:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
set -euo pipefail
CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:?}" && pwd)"
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
exec python3 "$MODEL_DIR/tools/run_e5f_joint_overnight_search.py" \
 --contract "${E5F_JOINT_CONTRACT:?}" --contract-sha256 "${E5F_JOINT_CONTRACT_SHA256:?}" \
 --mode "${E5F_JOINT_MODE:?}"
