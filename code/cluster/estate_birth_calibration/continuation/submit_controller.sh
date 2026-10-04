#!/usr/bin/env bash
# One controller only; it is queued behind successful completion of all ten parent tasks.
set -euo pipefail
remote=/scratch/td2248/projects/estate_birth_continuation_20261004_v1
ssh -o BatchMode=yes torch "cd '$remote' && test ! -e control/controller_submission.json && /share/apps/anaconda3/2025.06/bin/python verify_stage.py --host"
ssh -o BatchMode=yes torch "cd '$remote' && mkdir -p control && job=\$(sbatch --parsable --dependency=afterok:19127370 controller.sh) && job=\${job%%;*} && /share/apps/anaconda3/2025.06/bin/python - \"\$job\" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
receipt=dict(status='controller_queued_after_parent_success',controller_job_id=sys.argv[1],
  parent_array_job_id='19127370',dependency='afterok:19127370',
  stage_inventory_sha256=hashlib.sha256(Path('inventory.json').read_bytes()).hexdigest(),
  submitted_epoch=time.time(),no_auto_retry=True)
Path('control/controller_submission.json').write_text(json.dumps(receipt,indent=2)+chr(10))
print(json.dumps(receipt))
PY"
