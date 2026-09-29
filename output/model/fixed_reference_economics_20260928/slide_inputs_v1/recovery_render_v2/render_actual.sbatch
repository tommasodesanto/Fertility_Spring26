#!/usr/bin/env bash
#SBATCH --job-name=three-price-render-v2
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
root=/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v2
inputs=/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/solve_v1
module load anaconda3/2025.06
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
/share/apps/anaconda3/2025.06/bin/python "$root/source_v2/build_three_price_inputs.py" \
  --comparison "$inputs/comparison.csv" \
  --elasticities "$inputs/elasticities.csv" \
  --completed "$inputs/completed.json" \
  --output "$root/actual_output/recovery_18815133_v2"
