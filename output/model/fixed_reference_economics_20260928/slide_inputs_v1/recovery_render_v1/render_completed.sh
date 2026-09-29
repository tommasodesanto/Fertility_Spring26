#!/usr/bin/env bash
# Zero-solve rendering, submitted only after the approved recovery succeeds.
#SBATCH --job-name=three-price-render
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v1/render.log
#SBATCH --error=/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v1/render.err
set -euo pipefail
stage=/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v1
input=/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/solve_v1
printf '%s\n' "7e0d3a806dca7c0f27fce677921de5fef105a3aec3f11efa5e6beb9b5d06d766  $stage/source_v1/build_three_price_inputs.py" | sha256sum --check
module load anaconda3/2025.06
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
exec /share/apps/anaconda3/2025.06/bin/python "$stage/source_v1/build_three_price_inputs.py" \
  --comparison "$input/comparison.csv" --elasticities "$input/elasticities.csv" \
  --completed "$input/completed.json" --output "$stage/results/actual_v1"
