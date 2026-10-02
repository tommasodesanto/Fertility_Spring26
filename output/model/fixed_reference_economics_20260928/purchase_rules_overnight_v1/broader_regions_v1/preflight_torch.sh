#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/purchase_broader_regions_v1
mkdir -p "$remote/preflight"
for slot in 0 8; do
 SLURM_ARRAY_TASK_ID="$slot" REGIONS_PREFLIGHT=1 bash "$remote/launch_torch.sh"
 /share/apps/anaconda3/2025.06/bin/python - "$remote/preflight/slot_$slot" <<'PY'
import json,sys
from pathlib import Path
out=Path(sys.argv[1]);x=json.loads((out/'init/completed.json').read_text())
assert x['status']=='exact_initializer_passed_zero_lifecycle' and x['lifecycle_solves']==0
contract=json.loads((out/'init_region_contract.json').read_text())
assert contract['stage']=='init' and contract['slot'] in (0,8)
print(dict(slot=contract['slot'],arm=contract['arm'],status=x['status'],lifecycle_solves=0))
PY
done
