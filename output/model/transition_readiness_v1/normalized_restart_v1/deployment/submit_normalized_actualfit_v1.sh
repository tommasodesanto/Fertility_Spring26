#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/transition_readiness_v1/normalized_restart_v1
/share/apps/anaconda3/2025.06/bin/python - "$remote" <<'PY'
import hashlib,json,sys
from pathlib import Path
r=Path(sys.argv[1]);i=json.loads((r/'inventory.json').read_text())
assert hashlib.sha256((r/'inventory.json').read_bytes()).hexdigest()=='f3546476c38e339462321d5cb203950281de40a423ea2f8e443534840a885d1b'
for rel,h in i['files'].items():assert hashlib.sha256((r/'source'/rel).read_bytes()).hexdigest()==h,rel
assert hashlib.sha256((r/'floor_launch.sh').read_bytes()).hexdigest()==i['files']['code/model/experiments/transition_readiness/floor_launch.sh']
assert i['files']['code/model/experiments/transition_readiness/floor_runtime.py']=='cc51438a9a7bd0f92be3c1fe7681461515d21bb1ee44375ee7a51011489407af'
assert i['files']['output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json']=='5bbb824a1a8e2431efcfc227f65a5a2d2b8ce73a8b01001cffb74bf10d678a4e'
p=json.loads((r/'results/plan_preflight_v1/fit_preflight.json').read_text());assert p['status']=='PASS' and p['native_calls']==0
n=json.loads((r/'results/native_import_v1/preflight.json').read_text());assert n['policy_calls']==0 and n['native_setup_verified']
print('PASS frozen normalized inventory',len(i['files']))
PY
deadline=1790900687;remaining=$((deadline-$(date +%s)))
[[ "$remaining" -gt 120 ]] || exit 124
slurm_minutes=$(((remaining+59)/60))
sbatch --parsable --time="$slurm_minutes" --output="$remote/results/normalized_actualfit_v1_slurm.out" "$remote/floor_launch.sh" --mode fit --seconds "$remaining" --deadline-epoch "$deadline" --label normalized_actualfit_v1 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json
