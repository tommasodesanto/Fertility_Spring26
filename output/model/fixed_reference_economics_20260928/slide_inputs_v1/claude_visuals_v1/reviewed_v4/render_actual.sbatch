#!/usr/bin/env bash
#SBATCH --job-name=claude-visuals-v4
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:03:00
set -euo pipefail
root=/scratch/td2248/projects/fixed_reference_claude_visuals_20260929/reviewed_v4
module load anaconda3/2025.06
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
/share/apps/anaconda3/2025.06/bin/python "$root/render_claude_visuals_v4.py"
