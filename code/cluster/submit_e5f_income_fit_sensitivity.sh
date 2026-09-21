#!/usr/bin/env bash
# New immutable diagnostic panel; SUBMIT=1 submits smoke and afterok production.
set -euo pipefail
project_root="$(cd "$(dirname "$0")/../.." && pwd)"
local_root="$project_root/output/model/native_financing_diagnostic_20260919/specification_followup/quantification_v1/sensitivity"
remote_root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1
prior_root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1
mkdir -p "$local_root"
if [[ -e "$local_root/submission.json" || -e "$local_root/smoke_job.txt" ]]; then
  echo 'A submission receipt already exists; refusing duplicate launch.' >&2
  exit 2
fi
python3 - "$project_root" "$local_root" "$remote_root" "$prior_root" <<'PY'
from pathlib import Path
import json,sys,hashlib,datetime
project,local,remote,prior=map(Path,sys.argv[1:])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
driver=project/'code/model/tools/run_e5f_income_fit_sensitivity.py'
manifest={'status':'prepared_not_submitted','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'remote_root':str(remote),'driver_sha256':sha(driver),'prior_plan':str(prior/'plan.json'),'prior_plan_sha256':sha(project/'output/model/native_financing_diagnostic_20260919/overnight/plan.remote.json'),'selected_summary':str(prior/'production/summary.json'),'selected_summary_sha256':sha(project/'output/model/native_financing_diagnostic_20260919/overnight/final_search/summary.json'),'objective':'4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1','anchor_case':60,'anchor_loss':353.6588729140903,'full_probes':16,'half_probes':8,'anchor_evaluations':2,'selected_repetitions':2,'maximum_full_objective_evaluations':28,'maximum_stationary_solves_per_objective':8,'maximum_nested_stationary_solves':224,'observed_anchor_seconds':530.0652719750069,'expected_wall_minutes':[50,75],'workers':8,'case_timeout_seconds':900,'smoke_wall_seconds':3600,'production_wall_seconds':10800,'stop_submit_utc':'2026-09-21T06:00:00+00:00','finish_by_utc':'2026-09-21T11:30:00+00:00','interpretation':'controlled local fit sensitivity; no policy, global identification, convergence or adoption claim'}
(local/'launch_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for mode,wall in [('smoke',3300),('production',10200)]:
 extra='' if mode=='smoke' else f' --smoke-summary {remote}/smoke/smoke_summary.json'
 text=f'''#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
export NUMBA_CACHE_DIR={remote}/numba_cache_{mode}
export PYTHONPATH={prior}/code/model/tools
python - {remote} <<'CHECK'
from pathlib import Path
import hashlib,json,sys,datetime
r=Path(sys.argv[1]);m=json.loads((r/'launch_manifest.json').read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(r/'run_e5f_income_fit_sensitivity.py')==m['driver_sha256'],'new controller hash mismatch'
assert sha(m['prior_plan'])==m['prior_plan_sha256'],'parent plan hash mismatch'
assert sha(m['selected_summary'])==m['selected_summary_sha256'],'selected reference hash mismatch'
assert datetime.datetime.now(datetime.timezone.utc)<datetime.datetime.fromisoformat(m['finish_by_utc']),'absolute finish deadline already passed'
CHECK
remaining=$(python -c "import time,datetime;print(max(0,min({wall},int(datetime.datetime(2026,9,21,11,30,tzinfo=datetime.timezone.utc).timestamp()-time.time()))))")
[[ "$remaining" -gt 0 ]]
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "stage={mode} exit=%s\\n" "$rc" > {remote}/{mode}_failure.txt; fi' EXIT
timeout --signal=TERM --kill-after=20s "${{remaining}}s" python -B {remote}/run_e5f_income_fit_sensitivity.py --mode {mode} --plan {prior}/plan.json --selected-summary {prior}/production/summary.json --output {remote}/{mode} --workers 8 --case-timeout 900 --global-timeout "$remaining"{extra}
'''
 (local/f'{mode}.sbatch').write_text(text)
PY
bash -n "$local_root/smoke.sbatch"
bash -n "$local_root/production.sbatch"
if [[ "${DRY_RUN:-0}" == 1 ]]; then
  printf 'Prepared immutable panel at %s; no remote action.\n' "$local_root"
  exit 0
fi
# Refuse overwrite of any existing panel root.
ssh torch "mkdir '$remote_root'"
rsync -az "$project_root/code/model/tools/run_e5f_income_fit_sensitivity.py" "$local_root/launch_manifest.json" "$local_root/smoke.sbatch" "$local_root/production.sbatch" "torch:$remote_root/"
if [[ "${SUBMIT:-0}" != 1 ]]; then
  printf 'Staged only, not submitted: %s\n' "$remote_root"
  exit 0
fi
python3 - <<'PY'
import datetime
assert datetime.datetime.now(datetime.timezone.utc)<datetime.datetime(2026,9,21,6,tzinfo=datetime.timezone.utc),'submission cutoff passed'
PY
smoke=$(ssh torch "sbatch --parsable --account=torch_pr_570_general --job-name=income_fit_smoke --cpus-per-task=1 --mem=24G --time=01:00:00 --output=$remote_root/smoke_%j.out $remote_root/smoke.sbatch")
printf '%s\n' "$smoke" > "$local_root/smoke_job.txt"
production=$(ssh torch "sbatch --parsable --account=torch_pr_570_general --job-name=income_fit_panel --dependency=afterok:$smoke --kill-on-invalid-dep=yes --cpus-per-task=8 --mem=192G --time=03:00:00 --output=$remote_root/production_%j.out $remote_root/production.sbatch")
printf '%s\n' "$production" > "$local_root/production_job.txt"
python3 - "$local_root" "$remote_root" "$smoke" "$production" <<'PY'
from pathlib import Path
import sys,json,datetime
local,remote,smoke,production=sys.argv[1:]
r={'status':'submitted_production_gated_on_smoke','submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'remote_root':remote,'smoke_job':smoke,'production_job':production,'full_objective_evaluation_cap':28,'nested_stationary_solve_cap':224,'dependency':f'afterok:{smoke}','laptop_chaining_required':False,'finish_by_utc':'2026-09-21T11:30:00+00:00'}
Path(local,'submission.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
PY
