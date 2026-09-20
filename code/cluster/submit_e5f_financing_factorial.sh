#!/usr/bin/env bash
# Uses the immutable runtime staged by submit_e5f_income_overnight_search.sh.
set -euo pipefail
project_root="$(cd "$(dirname "$0")/../.." && pwd)"
local_root="$project_root/output/model/native_financing_diagnostic_20260919/overnight"
root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1
production_search=$(cat "$local_root/search_production_job.txt")
submit_matrix() {
 local mode="$1" dep="$2" name="$3" out="$4" family="$5"
 local dependency=""
 if [[ -n "$dep" ]]; then dependency="--dependency=afterok:$dep"; fi
 ssh torch sbatch --parsable --account=torch_pr_570_general ${dependency:+$dependency} --job-name="$name" --cpus-per-task=1 --mem=24G --time=03:00:00 --output="$root/${name}_%j.out" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
export NUMBA_CACHE_DIR="$root/numba_cache"
root="$root"
if [[ "$family" == refit ]]; then
 timeout --signal=TERM --kill-after=30s 10200s python -B "\$root/code/model/tools/run_e5f_financing_factorial.py" --mode "$mode" --plan "\$root/plan.json" --family refit_new_income --selection-summary "\$root/production/summary.json" --output "\$root/$out"
else
 timeout --signal=TERM --kill-after=30s 5000s python -B "\$root/code/model/tools/run_e5f_financing_factorial.py" --mode "$mode" --plan "\$root/plan.json" --family original --checkpoint /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz --expected-hash 3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993 --output "\$root/$out"
 timeout --signal=TERM --kill-after=30s 5000s python -B "\$root/code/model/tools/run_e5f_financing_factorial.py" --mode "$mode" --plan "\$root/plan.json" --family stationary_new_income --checkpoint /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_stationary_18040896/evaluation/raw/repetition_01/initial_state.pkl.gz --expected-hash 79f9dd5e56351bb6ced3a2b31d10657b1bfa8acd59c1c68436bfb671c99bb675 --output "\$root/$out"
fi
SBATCH
}
smoke=$(submit_matrix smoke '' finance_matrix_smoke factorial_smoke baseline)
printf '%s\n' "$smoke" > "$local_root/factorial_smoke_job.txt"
base=$(submit_matrix run "$smoke" finance_matrix_base factorial_baselines baseline)
printf '%s\n' "$base" > "$local_root/factorial_baselines_job.txt"
refit=$(submit_matrix run "$smoke:$production_search" finance_matrix_refit factorial_refit refit)
printf '%s\n' "$refit" > "$local_root/factorial_refit_job.txt"
printf 'matrix_smoke=%s matrix_baselines=%s matrix_refit=%s\n' "$smoke" "$base" "$refit"
