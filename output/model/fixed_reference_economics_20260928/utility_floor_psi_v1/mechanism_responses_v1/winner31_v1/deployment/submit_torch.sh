#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_winner31_responses_v1
cd "$remote"
[[ ! -e submission_receipt.json && ! -e results ]] || { echo 'Refusing duplicate submission'; exit 2; }
job=$(sbatch --parsable --time=00:20:00 launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY_SUBMIT'
import json,sys,time
from pathlib import Path
r=dict(status='submitted',job_id=sys.argv[1],submission_epoch=time.time(),deadline_basis='actual_launcher_start_plus_1200_seconds',cpu=1,memory_gib=24,threads=1,maximum_lifecycle_solves=6,case_seconds=300,total_seconds=1200,no_auto_retry=True,natural_support_certified=False,credit_ge_launched=False)
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
PY_SUBMIT
