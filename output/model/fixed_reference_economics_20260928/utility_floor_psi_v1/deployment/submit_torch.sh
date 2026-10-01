#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_psi_v2
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate submission'; exit 2; }
job=$(sbatch --parsable --array=0-7 launch_torch.sh)
job=${job%%;*}
python - "$job" <<'PY'
import json,sys,time
from pathlib import Path
r=dict(common_deadline_epoch=1790828667,status='submitted',array_job_id=sys.argv[1],chains=8,wall_seconds=7200,cpus_per_chain=1,memory_GiB_per_chain=24,maximum_objective_calls_per_chain=200,free_psi=True,profiles='original_baseline_weights',experimental_not_adopted=True,submission_epoch=time.time(),no_auto_retry=True)
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
PY
