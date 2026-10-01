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
assert i['files']['code/model/experiments/transition_readiness/one_shock_floor.py']=='a17f14da40945a394e75cbfee4f46009e18fe7e28d824f2d23b9cdbdd83296b2'
print('PASS remote frozen inventory',len(i['files']))
PY
sbatch --parsable --time=00:11:00 --output="$remote/results/observer_debug_v1_slurm.out" "$remote/floor_launch.sh" --mode debug --seconds 600 --label observer_debug_v1 --receipt /work/transition_runs/smoke_v3/native_reference/repeat_0/native_solve_unverified.json --receipt-sha256 ab223ec92324949c795bb8f38f3680a82d59a2ed3a465733a130a25682ec95e6
