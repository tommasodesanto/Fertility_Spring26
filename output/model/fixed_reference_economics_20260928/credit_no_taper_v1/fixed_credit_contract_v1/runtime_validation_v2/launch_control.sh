#!/usr/bin/env bash
#SBATCH --job-name=e5f_credit_control_v2
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
set -euo pipefail
exec "$(dirname "$0")/launch_runtime_validation.sh" control
