#!/usr/bin/env bash
# Invoke only after explicit root release; paired smoke gate mandatory.
set -euo pipefail
remote=/scratch/td2248/projects/estate_birth_calibration_20261003_v3
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate production submission'; exit 2; }
/share/apps/anaconda3/2025.06/bin/python verify_stage.py --host
/share/apps/anaconda3/2025.06/bin/python verify_smoke_gate.py "$remote"
/share/apps/anaconda3/2025.06/bin/python - <<'PY'
import os
s=os.statvfs('.')
assert s.f_bavail*s.f_frsize>=5200*1024**3,'Insufficient planned shared scratch reserve (5200GiB)'
PY
job=$(sbatch --parsable --array=0-9%10 --time=06:00:00 --export=ALL,ESTATE_RUN_MODE=production launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
receipt=dict(status='production_submitted',job_id=sys.argv[1],array='0-9%10',chains=10,account='torch_pr_570_general',
cpus_per_chain=1,memory_GiB_per_chain=24,wall_seconds_per_chain=21600,maximum_objective_calls_per_chain=500,
final_native_reserve_seconds=1800,no_auto_retry=True,submission_epoch=time.time(),
inventory_sha256=hashlib.sha256(Path('inventory.json').read_bytes()).hexdigest())
Path('submission_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
PY
