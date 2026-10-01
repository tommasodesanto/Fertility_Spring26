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
assert i['files']['code/model/experiments/transition_readiness/one_shock_floor.py']=='658859f10ec709d92f07bc599d3d36144c672827a99ffdd16bc4cee6892affef'
assert i['files']['output/model/transition_readiness_v1/current_floor_preparation/deployment/smoke_plan.json']=='7e046a013d94b68eb23a3ddb86572c79b146d8f1a1c374045b281bdd86154294'
print('PASS remote frozen inventory',len(i['files']))
PY
sbatch --parsable --time=01:31:00 --output="$remote/results/smoke_v4_slurm.out" "$remote/floor_launch.sh" --mode smoke --seconds 5400 --label smoke_v4 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_floor_preparation/deployment/smoke_plan.json
