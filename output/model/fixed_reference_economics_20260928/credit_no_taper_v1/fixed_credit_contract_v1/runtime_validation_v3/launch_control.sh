#!/usr/bin/env bash
#SBATCH --job-name=e5f_credit_control_v3
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
set -euo pipefail
exec "/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v3/source/launch_runtime_validation.sh" control
