#!/usr/bin/env bash
# Stage a separate, hash-pinned coordinate poll; never overwrite the pilot runtime.
set -euo pipefail
project_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_search_v1
local_root="$project_root/output/model/native_financing_diagnostic_20260919/earnings_candidate"
python3 - "$project_root" "$run_root" <<'PY'
import hashlib,json,sys
from pathlib import Path
root,remote=map(Path,sys.argv[1:])
local=root/'output/model/native_financing_diagnostic_20260919/earnings_candidate'
p=json.loads((local/'calibration_plan.remote.json').read_text())
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for key,name in [('constructor','build_persistent_transitory_income_candidate.py'),('adapter','run_e5f_income_candidate_calibration.py'),('controller','run_e5f_income_candidate_search.py')]:
    p[key+'_path']=str(remote/'code/model/tools'/name)
    p[key+'_sha256']=sha(root/'code/model/tools'/name)
p['candidate_json']=str(remote/'earnings_candidate/candidate.json')
p['candidate_json_sha256']=sha(local/'candidate.json')
p['plan_path']=str(remote/'earnings_candidate/search_plan.json')
score_path=root/'output/model/native_financing_diagnostic_20260919/income_stationary_18040896/evaluation/scored_repetition_01/score.json'
s=json.loads(score_path.read_text())
p['incumbent_score_sha256']=sha(score_path)
sig=[(r['restriction_id'],r['target'],r.get('actual_weight'),r.get('scored')) for r in s['target_fit']]
p['target_signature']=hashlib.sha256(json.dumps(sig,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
p['parameter_bounds']={r['parameter']:[r['lower'],.99 if r['parameter']=='beta_annual' else r['upper']] for r in s['parameters'] if r.get('structural_coordinate')}
p['incumbent_score']=str(remote.parent/'income_stationary_18040896/evaluation/scored_repetition_01/score.json')
p['incumbent_checkpoint']=str(remote.parent/'income_stationary_18040896/evaluation/raw/repetition_01/initial_state.pkl.gz')
p['incumbent_checkpoint_sha256']='79f9dd5e56351bb6ced3a2b31d10657b1bfa8acd59c1c68436bfb671c99bb675'
p.update(search_max_proposals=16,search_workers=4,final_repetitions=2,
         search_seconds=2100,case_timeout=900,verification_seconds=1200,
         hard_timeout_seconds=3420,slurm_seconds=3600,
         estimated_trial_seconds=[200,500],estimated_native_solves=[36,144],
         stop_condition='One finite coordinate poll; stop dispatch with under 500s remaining; repeat selected point twice; abort on source/income/target mismatch.',
         interpretation='bounded local diagnostic refit, not convergence or adoption',
         status='authorized_bounded_income_refit')
(local/'search_plan.remote.json').write_text(json.dumps(p,indent=2,sort_keys=True)+'\n')
PY
ssh torch "mkdir -p '$run_root/code/model/tools' '$run_root/earnings_candidate'"
rsync -az "$project_root/code/model/tools/run_e5f_income_candidate_search.py" "$project_root/code/model/tools/run_e5f_income_candidate_calibration.py" "$project_root/code/model/tools/build_persistent_transitory_income_candidate.py" "torch:$run_root/code/model/tools/"
rsync -az "$local_root/candidate.json" "torch:$run_root/earnings_candidate/"
rsync -az "$local_root/search_plan.remote.json" "torch:$run_root/earnings_candidate/search_plan.json"
if [[ "${SUBMIT:-0}" != 1 ]]; then
    printf 'Staged %s; set SUBMIT=1 to launch.\n' "$run_root"
    exit 0
fi
ssh torch "sbatch --parsable --account=torch_pr_570_general --job-name=e5f_income_search --cpus-per-task=4 --mem=96G --time=01:00:00 --output=$run_root/search_%j.out" <<'SBATCH' | tee "$local_root/search_submission.txt"
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1 NUMBA_DISABLE_JIT=0
run_root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_search_v1
export NUMBA_CACHE_DIR="$run_root/numba_cache"
plan="$run_root/earnings_candidate/search_plan.json"
controller="$run_root/code/model/tools/run_e5f_income_candidate_search.py"
output="$run_root/results_${SLURM_JOB_ID}"
incumbent=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_stationary_18040896/evaluation/scored_repetition_01/score.json
cd "$run_root"
# Zero-solve actual-checkpoint income audit and two-repetition wrapper contract.
python -B "$controller" --mode preflight --plan "$plan" --output "${output}_preflight"
set +e
timeout --signal=TERM --kill-after=30s 3420s python -B "$controller" --mode search --plan "$plan" --output "$output" --incumbent-score "$incumbent"
result=$?
printf '{"job_id":"%s","exit_code":%d,"scope":"at most 16 coordinate proposals plus two native repetitions of selected point; not a convergence claim"}\n' "$SLURM_JOB_ID" "$result" > "${output}_job_receipt.json"
exit "$result"
SBATCH
