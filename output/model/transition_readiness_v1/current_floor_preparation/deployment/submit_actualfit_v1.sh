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
assert i['files']['code/model/experiments/transition_readiness/one_shock_floor.py']=='880fe888b3863e45fa10ac59538868f02033b8289a1f5ba6d7c500d0ffae203c'
assert i['files']['output/model/transition_readiness_v1/current_floor_preparation/deployment/fit_plan.json']=='4a83a2b5496362876a30c68dcacd485c6f30dcde95aaaf76c4c70c0edb89dd53'
pre=json.loads((r/'results/actualfit_preflight_v1/fit_preflight.json').read_text())
assert pre['status']=='exact_mounted_preparation_preflight_passed' and pre['native_calls']==0
print('PASS remote frozen inventory',len(i['files']))
PY
sbatch --parsable --time=05:01:00 --output="$remote/results/actualfit_v1_slurm.out" "$remote/floor_launch.sh" --mode fit --seconds 18000 --label actualfit_v1 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_floor_preparation/deployment/fit_plan.json
