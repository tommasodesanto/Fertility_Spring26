#!/bin/bash
#SBATCH --job-name=e5f_price_jacobian
#SBATCH --output=logs/price_jacobian_%j.out
#SBATCH --error=logs/price_jacobian_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a
printf '%s  %s\n' 'a002192b0a3cb9bbc561ca7d7b12185a10e1a094cbc55ff860d18f0123ad6d5e' 'code/model/tools/collect_e5f_matched_pf_price_jacobian.py' | sha256sum --check --status
python code/model/tools/collect_e5f_matched_pf_price_jacobian.py --anchor /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/path_anchor_02/sequential --probes /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/price_probe_01/sequential/{0..11} --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a/output/jacobian_sequential.json
