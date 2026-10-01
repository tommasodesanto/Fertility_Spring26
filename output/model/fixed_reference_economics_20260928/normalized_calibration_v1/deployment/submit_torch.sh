#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/normalized_calibration_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate submission'; exit 2; }
smoke=$(sbatch --parsable --array=0-0%1 --time=02:00:00 launch_smoke_torch.sh)
smoke=${smoke%%;*}
prod=$(sbatch --parsable --dependency="afterok:$smoke" --array=0-23%24 --time=02:00:00 launch_torch.sh)
prod=${prod%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$smoke" "$prod" <<'PY'
import json,sys,time
from pathlib import Path
r=dict(status='submitted_dependency_gated',smoke_job_id=sys.argv[1],array_job_id=sys.argv[2],chains=24,dependency=f'afterok:{sys.argv[1]}',account='torch_pr_570_general',cpus_per_chain=1,memory_GiB_per_chain=24,threads_per_chain=1,actual_start_seconds_per_chain=7200,maximum_objective_calls_per_chain=100,final_native_reserve_seconds=900,objective_weights='original_base_weights_only',experimental_not_adopted=True,no_auto_retry=True,submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
PY
