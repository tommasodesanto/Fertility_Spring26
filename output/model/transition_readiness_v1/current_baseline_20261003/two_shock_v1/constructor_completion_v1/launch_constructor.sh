#!/usr/bin/env bash
#SBATCH --job-name=constructor_completion
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cl
set -euo pipefail
module load anaconda3/2025.06
exec /share/apps/anaconda3/2025.06/bin/python /scratch/td2248/projects/current_estate_two_shock_20261004_v4_constructor_completion_v1/prepare_two_shock_torch.py run-constructor-completion --packet /scratch/td2248/projects/current_estate_two_shock_20261004_v4_constructor_completion_v1 --inventory-sha256 fe7d5ea89436901a36c233f821931ec4ce75f56b4dcd26ae69e74345d05aba1a
