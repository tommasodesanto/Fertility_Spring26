#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/transition_readiness_v1/normalized_resumed_fit_v1
/share/apps/anaconda3/2025.06/bin/python - "$remote" <<'PY'
import json,hashlib,sys
from pathlib import Path
r=Path(sys.argv[1]);inv=json.loads((r/'inventory.json').read_text());assert hashlib.sha256((r/'inventory.json').read_bytes()).hexdigest()=='bedfc49049521f0167569eecc706606eccc5e6938504e7c1bb21e26e9f78e62a'
for rel,h in inv['files'].items():assert hashlib.sha256((r/'source'/rel).read_bytes()).hexdigest()==h,rel
assert hashlib.sha256((r/'floor_launch.sh').read_bytes()).hexdigest()==inv['files']['code/model/experiments/transition_readiness/floor_launch.sh']
p=r/'source/output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/resume_plan.json';assert hashlib.sha256(p.read_bytes()).hexdigest()=='f33912119787da38de6c72cf238b5c6366bac90763f946494b37f6035910a05f';plan=json.loads(p.read_text())
assert plan['source_files']['controller']['sha256']=='678b602da22df0fb9a39cd8cf4772769ee4b6f06f8432976099d24712674c051'
n=json.loads((r/'results/normalized_resume_native_import_v1/preflight.json').read_text());assert n['policy_calls']==0 and n['native_setup_verified'] and n['identity']==plan['identity']
a=json.loads((r/'results/normalized_resume_plan_preflight_v1/fit_preflight.json').read_text());assert a['status']=='exact_mounted_preparation_preflight_passed' and a['native_calls']==0
print('PASS frozen resumed inventory',len(inv['files']))
PY
# One atomic reservation; an unknown submission outcome never triggers a duplicate.
mkdir "$remote/submission.lock"
job=$(sbatch --parsable --time=02:01:00 --output="$remote/results/normalized_resumed_fit_v1_slurm.out" "$remote/floor_launch.sh" --mode fit --seconds 7200 --label normalized_resumed_fit_v1 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/normalized_restart_v1/resume_preparation/deployment/resume_plan.json)
/share/apps/anaconda3/2025.06/bin/python - "$remote" "$job" <<'PY'
import json,os,sys,time
from pathlib import Path
p=Path(sys.argv[1])/'submission_receipt.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(job_id=sys.argv[2],submitted_epoch=time.time(),sole_submitter='lead',label='normalized_resumed_fit_v1',numerical_seconds=7200,slurm_seconds=7260,maximum_policy_calls=1214),indent=2)+'\n');os.replace(t,p)
PY
printf '%s\n' "$job"
