#!/bin/bash
#SBATCH --account=torch_pr_570_general
#SBATCH --job-name=e5f_rebated_joint_search
#SBATCH --cpus-per-task=6
#SBATCH --mem=72G
#SBATCH --time=03:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/search_joint_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/search_joint_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/output/cache/numba
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_source:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/initial_v2:/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/md_exact_loop/inputs:/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/code/model/tools:/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/code/model
python -B -m unittest test_run_e5f_rebated_initial_search -q
python -B /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_source/run_e5f_rebated_initial_search.py run --helper /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/initial_v2/run_e5f_rebated_initial_overnight.py --helper-sha256 d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790 --joint /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_source/run_e5f_joint_rebated_initial_probe.py --joint-sha256 9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb --scored /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_source/run_e5f_joint_rebated_initial_scored.py --scored-sha256 9da35b55466d74dc10a6a85ff91dabe887a807aab68417eabf63a9d7a42ca967 --template /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/capped_beta_099_20260911 --smoke /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_initial_smoke --output /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/initial_search_joint --workers 6 --wall-seconds 10800
