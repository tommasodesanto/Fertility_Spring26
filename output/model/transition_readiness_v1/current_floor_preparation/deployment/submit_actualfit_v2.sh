#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/transition_readiness_v1/current_floor
/share/apps/anaconda3/2025.06/bin/python - "$remote" <<'PY'
import hashlib,json,sys
from pathlib import Path
r=Path(sys.argv[1]);i=json.loads((r/'inventory.json').read_text())
for rel,digest in i['files'].items():assert hashlib.sha256((r/'source'/rel).read_bytes()).hexdigest()==digest,rel
assert hashlib.sha256((r/'floor_launch.sh').read_bytes()).hexdigest()==i['files']['code/model/experiments/transition_readiness/floor_launch.sh']
assert i['files']['code/model/experiments/transition_readiness/floor_runtime.py']=='f5ad7bf6c6d86962b9bc436b3e2d82951d7d7d3aac699b69abf340b0c0a508fd'
assert i['files']['code/model/experiments/transition_readiness/one_shock_floor.py']=='d011b35432706ceb011e605e3a4eaf7dcbb9cf95e76b5f3faa2ffb6875e7b8fb'
assert i['files']['output/model/transition_readiness_v1/current_floor_preparation/deployment/fit_v2_plan.json']=='6b6b28afbf798562f2f28fd8c38eb33d844a36363bc53c45d4b35b87aeb34945'
pre=json.loads((r/'results/actualfit_preflight_v2/fit_preflight.json').read_text())
assert pre['status']=='exact_mounted_preparation_preflight_passed' and pre['native_calls']==0
print('PASS remote frozen inventory',len(i['files']))
PY
deadline=1790900687
remaining=$((deadline-$(date +%s)))
[[ "$remaining" -gt 120 ]] || { echo 'Original deadline exhausted'; exit 124; }
slurm_minutes=$(((remaining+59)/60+1))
sbatch --parsable --time="$slurm_minutes" --output="$remote/results/actualfit_v2_slurm.out" "$remote/floor_launch.sh" --mode fit --seconds "$remaining" --deadline-epoch "$deadline" --label actualfit_v2 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_floor_preparation/deployment/fit_v2_plan.json
