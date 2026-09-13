#!/bin/bash
#SBATCH --account=torch_pr_570_general
#SBATCH --job-name=e5f_cache_probe
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/cache_probe_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/cache_probe_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/cache_source_v2:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/history_source:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/code/model/tools:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/code/model
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/output/cache/numba
python -B -m unittest test_e5f_exact_policy_cache -q
python -B /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/cache_source_v2/run_e5f_policy_cache_probe.py --manifest /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/history_manifest_joint.json --output /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/policy_cache_probe_v2 --cache-sha256 d51bbd13026288db0194b44b420bb49ff6970588f40a906242978d49038a5d6f
