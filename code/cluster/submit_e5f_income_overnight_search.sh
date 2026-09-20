#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1
local_root="$project_root/output/model/native_financing_diagnostic_20260919/overnight"
mkdir -p "$local_root"
python3 - "$project_root" "$run_root" "$local_root" <<'PY'
import hashlib,json,sys
from pathlib import Path
root,remote,local=map(Path,sys.argv[1:]);base=root/'output/model/native_financing_diagnostic_20260919/earnings_candidate'
p=json.loads((base/'search_plan.remote.json').read_text());sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for key,name in [('constructor','build_persistent_transitory_income_candidate.py'),('adapter','run_e5f_income_candidate_calibration.py'),('controller','run_e5f_income_candidate_search.py'),('overnight_controller','run_e5f_income_overnight_search.py'),('factorial_controller','run_e5f_financing_factorial.py')]:
 p[key+'_path']=str(remote/'code/model/tools'/name);p[key+'_sha256']=sha(root/'code/model/tools'/name)
p['factorial_helper_sha256']={n:sha(root/'code/model/tools'/n) for n in ['run_e5f_native_financing_diagnostic.py','run_e5f_native_rental_access_diagnostic.py','run_e5f_native_income_cohort_diagnostic.py','build_e5f_native_financing_report.py']}
p['candidate_json']=str(remote/'earnings_candidate/candidate.json');p['plan_path']=str(remote/'plan.json')
p.update(status='authorized_overnight_diagnostic',overnight_max_proposals=96,overnight_workers=8,overnight_search_seconds=16200,overnight_case_timeout=900,overnight_verification_seconds=1800,overnight_slurm_seconds=21600,estimated_point_seconds=400,estimated_native_solves_max=784,estimated_search_wall_minutes=110,search_scope='finite joint near/wide proposals; no convergence/adoption claim',factorial_scope='three binary finance/rental margins within each checkpoint family; fixed prices and preferences; no GE claim')
(local/'plan.remote.json').write_text(json.dumps(p,indent=2,sort_keys=True)+'\n')
PY
ssh torch "mkdir -p '$run_root/code/model/tools' '$run_root/earnings_candidate'"
files=(build_persistent_transitory_income_candidate.py run_e5f_income_candidate_calibration.py run_e5f_income_candidate_search.py run_e5f_income_overnight_search.py run_e5f_financing_factorial.py run_e5f_native_financing_diagnostic.py run_e5f_native_rental_access_diagnostic.py run_e5f_native_income_cohort_diagnostic.py build_e5f_native_financing_report.py)
for file in "${files[@]}"; do rsync -az "$project_root/code/model/tools/$file" "torch:$run_root/code/model/tools/"; done
rsync -az "$project_root/output/model/native_financing_diagnostic_20260919/earnings_candidate/candidate.json" "torch:$run_root/earnings_candidate/"
rsync -az "$local_root/plan.remote.json" "torch:$run_root/plan.json"
if [[ "${SUBMIT:-0}" != 1 ]]; then printf 'Staged %s\n' "$run_root"; exit 0; fi
smoke_job=$(ssh torch "sbatch --parsable --account=torch_pr_570_general --dependency=afterok:18047156 --job-name=income_joint_smoke --cpus-per-task=2 --mem=48G --time=01:00:00 --output=$run_root/smoke_%j.out" <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1
export NUMBA_CACHE_DIR="$root/numba_cache"
prior=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_search_v1/results_18047156/summary.json
timeout --signal=TERM --kill-after=30s 3420s python -B "$root/code/model/tools/run_e5f_income_overnight_search.py" --mode smoke --plan "$root/plan.json" --prior-summary "$prior" --output "$root/smoke"
SBATCH
)
printf '%s\n' "$smoke_job" > "$local_root/search_smoke_job.txt"
production_job=$(ssh torch "sbatch --parsable --account=torch_pr_570_general --dependency=afterok:$smoke_job --job-name=income_joint_night --cpus-per-task=8 --mem=192G --time=06:00:00 --output=$run_root/production_%j.out" <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1
export NUMBA_CACHE_DIR="$root/numba_cache"
prior=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_search_v1/results_18047156/summary.json
timeout --signal=TERM --kill-after=30s 21000s python -B "$root/code/model/tools/run_e5f_income_overnight_search.py" --mode search --plan "$root/plan.json" --prior-summary "$prior" --smoke-summary "$root/smoke/summary.json" --output "$root/production"
SBATCH
)
printf '%s\n' "$production_job" > "$local_root/search_production_job.txt"
printf 'smoke=%s production=%s\n' "$smoke_job" "$production_job"
