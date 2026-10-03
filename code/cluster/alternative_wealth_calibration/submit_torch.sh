#!/usr/bin/env bash
# Production submission is a separate lead-reviewed action, after native smoke.
set -euo pipefail
remote=/scratch/td2248/projects/alternative_wealth_calibration_20261003_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate production submission'; exit 2; }
/share/apps/anaconda3/2025.06/bin/python verify_stage.py --host
/share/apps/anaconda3/2025.06/bin/python verify_smoke_gate.py "$remote"
job=$(sbatch --parsable --array=0-9%10 --time=06:00:00 --export=ALL,WEALTH_RUN_MODE=production launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
receipt=dict(status='production_submitted',job_id=sys.argv[1],array='0-9%10',chains=10,
             account='torch_pr_570_general',cpus_per_chain=1,memory_GiB_per_chain=24,
             wall_seconds_per_chain=21600,maximum_objective_calls_per_chain=500,
             final_native_reserve_seconds=1800,no_auto_retry=True,
             inventory_sha256=hashlib.sha256(Path('inventory.json').read_bytes()).hexdigest(),
             submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
PY
