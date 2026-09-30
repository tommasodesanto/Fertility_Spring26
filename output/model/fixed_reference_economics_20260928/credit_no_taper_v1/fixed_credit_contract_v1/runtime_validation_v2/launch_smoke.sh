#!/usr/bin/env bash
#SBATCH --job-name=e5f_credit_smoke_v2
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
set -euo pipefail
exec "$(dirname "$0")/launch_runtime_validation.sh" smoke
