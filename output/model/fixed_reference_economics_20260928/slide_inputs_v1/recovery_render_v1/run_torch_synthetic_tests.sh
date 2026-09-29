#!/usr/bin/env bash
# Synthetic-only schema and chart checks; no economic result is rendered.
#SBATCH --job-name=three-price-test
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail

stage=/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v1
module load anaconda3/2025.06
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1

/share/apps/anaconda3/2025.06/bin/python "$stage/source_v1/build_three_price_inputs.py" --help > "$stage/results/cli_help.txt"
/share/apps/anaconda3/2025.06/bin/python "$stage/source_v1/test_three_price_inputs.py" "$stage/results/synthetic_test_only"
