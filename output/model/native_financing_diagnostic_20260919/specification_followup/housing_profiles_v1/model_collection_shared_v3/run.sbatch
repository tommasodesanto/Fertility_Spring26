#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --job-name=housing_saved_v3
#SBATCH --chdir=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/housing_profiles_shared_v3
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/housing_profiles_shared_v3/slurm_%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONUNBUFFERED=1
python - <<'VERIFY'
import hashlib,json,pathlib
m=json.load(open('launch_manifest.json'))
for name,pin in m['inputs'].items():
 assert hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest()==pin,name
VERIFY
timeout 80 python -B collect_e5f_saved_housing_profiles.py --checkpoint /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz --source-root /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2 --summary original_summary.json --family original --output original.json
timeout 80 python -B collect_e5f_saved_housing_profiles.py --checkpoint /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_stationary_18040896/evaluation/raw/repetition_01/initial_state.pkl.gz --source-root /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2 --summary stationary_new_income_summary.json --family stationary_new_income --output stationary_new_income.json
timeout 80 python -B collect_e5f_saved_housing_profiles.py --checkpoint /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/production/selected_verification/evaluation/raw/repetition_02/initial_state.pkl.gz --source-root /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2 --summary refit_new_income_summary.json --family refit_new_income --output refit_new_income.json
