#!/usr/bin/env bash
# Explicit production submission after lead review; no dependency auto-submit.
set -euo pipefail
remote=/scratch/td2248/projects/soft_timing_calibration_20261002_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate production submission'; exit 2; }
/share/apps/anaconda3/2025.06/bin/python verify_stage.py --host
/share/apps/anaconda3/2025.06/bin/python verify_smoke_gate.py "$remote"
job=$(sbatch --parsable --array=0-7%8 --time=06:00:00 --export=ALL,SOFT_RUN_MODE=production launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY'
import json,sys,time
from pathlib import Path
receipt=dict(status='production_submitted',job_id=sys.argv[1],array='0-7%8',
             arms={'original':'0-3','alternative':'4-7'},chains_per_arm=4,
             account='torch_pr_570_general',cpus_per_chain=1,memory_GiB_per_chain=24,
             wall_seconds_per_chain=21600,maximum_objective_calls_per_chain=250,
             final_native_reserve_seconds=1800,no_auto_retry=True,
             submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
PY
