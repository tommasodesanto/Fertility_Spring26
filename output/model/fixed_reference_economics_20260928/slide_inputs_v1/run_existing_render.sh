#!/usr/bin/env bash
# Zero-solve renderer for already-passed borrowing and supply results.
#SBATCH --job-name=fixed-econ-fig
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail

stage=/scratch/td2248/projects/fixed_reference_economic_figures_20260929
credit=/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1
ge=/scratch/td2248/projects/fixed_reference_credit_ge_20260929/ge_results/solve_v1/selected_repeat
supply=/scratch/td2248/projects/fixed_reference_supply_20260929/results_v2/replay
module load anaconda3/2025.06
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1

/share/apps/anaconda3/2025.06/bin/python "$stage/source_v1/build_slide_inputs.py" \
  --reference-receipt "$credit/control/receipt.json" \
  --grid-receipt "$credit/grid_control/receipt.json" \
  --credit-receipt "$credit/credit/receipt.json" \
  --ge-receipt "$ge/receipt.json" \
  --reference-fit "$credit/control/target_fit.csv" \
  --grid-fit "$credit/grid_control/target_fit.csv" \
  --credit-fit "$credit/credit/target_fit.csv" \
  --ge-fit "$ge/target_fit.csv" \
  --supply-comparison "$supply/comparison.csv" \
  --supply-verification "$supply/verification.json" \
  --output "$stage/results/output_v1"
