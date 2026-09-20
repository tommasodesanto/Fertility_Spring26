#!/usr/bin/env bash
# Isolated 48-cell dose experiment; smoke-gated and checkpointed on Torch.
set -euo pipefail
project_root="$(cd "$(dirname "$0")/../.." && pwd)"
family="${FAMILY:-original}"
tag="${DOSE_TAG:-finance_dose_v1}"
population_source="${POPULATION_SOURCE:-stationary}"
case "$family" in original|stationary_new_income|refit_new_income) ;; *) printf 'Invalid family\n' >&2; exit 2 ;; esac
case "$population_source" in stationary|saved_evaluation) ;; *) printf 'Invalid population source\n' >&2; exit 2 ;; esac
[[ "$tag" =~ ^[a-zA-Z0-9_-]+$ ]] || { printf 'Invalid DOSE_TAG\n' >&2; exit 2; }
root="/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/$tag"
local_root="$project_root/output/model/native_financing_diagnostic_20260919/overnight/$tag"
# Never overwrite a snapshot that may already have queued or running jobs.
ssh torch "mkdir '$root'" || { printf 'Use a fresh DOSE_TAG; snapshot exists or creation failed.\n' >&2; exit 2; }
mkdir -p "$local_root"
python3 - "$project_root" "$root" "$local_root" "$population_source" <<'PY'
import json,hashlib,sys
from pathlib import Path
repo,remote,local=map(Path,sys.argv[1:4]);base=repo/'output/model/native_financing_diagnostic_20260919/overnight/plan.remote.json'
p=json.loads(base.read_text());sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for k in ['constructor','adapter','controller','overnight_controller','factorial_controller']:
 name=Path(p[k+'_path']).name;p[k+'_path']=str(remote/'code/model/tools'/name);p[k+'_sha256']=sha(repo/'code/model/tools'/name)
p['factorial_helper_sha256']={name:sha(repo/'code/model/tools'/name) for name in p['factorial_helper_sha256']}
p['candidate_json']=str(remote/'earnings_candidate/candidate.json');p['plan_path']=str(remote/'plan.json')
p['dose_design']={'phi':[.8,.9,.95,1.],'lambda_four_year_income':[0.,.25,1.,5.],'rental_cap':[6.,8.,10.],'cells':48,'extra_controls':2,'case_seconds':600,'total_seconds':18000,'slurm_seconds':21600,'threads':1,'estimated_case_seconds':[60,120],'estimated_wall_minutes':[50,100]}
p['dose_design']['population_source']=sys.argv[4]
(local/'plan.remote.json').write_text(json.dumps(p,indent=2,sort_keys=True)+'\n')
PY
ssh torch "mkdir -p '$root/code/model/tools' '$root/earnings_candidate'"
files=(build_persistent_transitory_income_candidate.py run_e5f_income_candidate_calibration.py run_e5f_income_candidate_search.py run_e5f_income_overnight_search.py run_e5f_financing_factorial.py run_e5f_native_financing_diagnostic.py run_e5f_native_rental_access_diagnostic.py run_e5f_native_income_cohort_diagnostic.py build_e5f_native_financing_report.py)
paths=(); for file in "${files[@]}"; do paths+=("$project_root/code/model/tools/$file"); done
rsync -az "${paths[@]}" "torch:$root/code/model/tools/"
rsync -az "$project_root/output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json" "torch:$root/earnings_candidate/"
rsync -az "$local_root/plan.remote.json" "torch:$root/plan.json"
if [[ "${SUBMIT:-0}" != 1 ]]; then printf 'Staged %s\n' "$root"; exit 0; fi
submit() {
 local mode="$1" dep="$2" output="$3"
 local wall="06:00:00"; if [[ "$mode" == smoke ]]; then wall="01:00:00"; fi
 local depflag=""; if [[ -n "$dep" ]]; then depflag="--dependency=afterok:$dep"; fi
 ssh torch sbatch --parsable --account=torch_pr_570_general ${depflag:+$depflag} --job-name="finance_dose_${family}_${mode}" --cpus-per-task=1 --mem=24G --time="$wall" --output="$root/${family}_${mode}_%j.out" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
export NUMBA_CACHE_DIR="$root/numba_cache"
args=(--mode "$mode" --design dose --total-seconds 18000 --population-source "$population_source" --plan "$root/plan.json" --family "$family" --output "$root/$output")
if [[ "$family" == original ]]; then
 args+=(--checkpoint /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz --expected-hash 3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993)
elif [[ "$family" == stationary_new_income ]]; then
 args+=(--checkpoint /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_stationary_18040896/evaluation/raw/repetition_01/initial_state.pkl.gz --expected-hash 79f9dd5e56351bb6ced3a2b31d10657b1bfa8acd59c1c68436bfb671c99bb675)
else
 args+=(--selection-summary /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/production/summary.json)
fi
timeout --signal=TERM --kill-after=30s 21000s python -B "$root/code/model/tools/run_e5f_financing_factorial.py" "\${args[@]}"
SBATCH
}
smoke=$(submit smoke "${AFTER_JOB:-}" "${family}_smoke")
printf '%s\n' "$smoke" > "$local_root/${family}_smoke_job.txt"
run=$(submit run "$smoke" "${family}_production")
printf '%s\n' "$run" > "$local_root/${family}_production_job.txt"
printf 'family=%s smoke=%s production=%s\n' "$family" "$smoke" "$run"
