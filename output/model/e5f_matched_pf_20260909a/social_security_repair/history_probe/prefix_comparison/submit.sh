#!/bin/bash
#SBATCH --job-name=e5f_prefix_compare
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_paygo_prefix_527ab397/prefix_comparison_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_paygo_prefix_527ab397/prefix_comparison_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_paygo_prefix_527ab397
export PYTHONPATH="$PWD/code/model/tools:$PWD/code/model"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/e5f_prefix_mpl_${SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"
python prefix_comparison/run_prefix_comparison.py --contract prefix_comparison/contract.json --contract-sha256 698b77087e08964c04d4331dad5c99f9fdec063d04624b5d9ee3ea37a380ddc1 --source-root "$PWD" --output "output/prefix_comparison_${SLURM_JOB_ID}"
