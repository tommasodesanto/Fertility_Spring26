#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/purchase_rules_overnight_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate submission'; exit 2; }
smoke=$(sbatch --parsable --array=0,24 --time=00:20:00 --export=ALL,CALIBRATION_SMOKE=1 launch_torch.sh)
smoke=${smoke%%;*}
production=$(sbatch --parsable --dependency="afterok:$smoke" --array=0-47%48 --time=04:00:00 --export=ALL,CALIBRATION_SMOKE=0 launch_torch.sh)
production=${production%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$smoke" "$production" <<'PY'
import json,sys,time
from pathlib import Path
receipt=dict(status='submitted_smoke_dependency_gated',smoke_array_job_id=sys.argv[1],production_array_job_id=sys.argv[2],arms={'hard':'0-23','quarter':'24-47'},production_chains=48,dependency='afterok:'+sys.argv[1],cpus_per_chain=1,memory_GiB_per_chain=24,threads_per_chain=1,wall_seconds_per_chain=14400,maximum_objective_calls_per_chain=250,final_native_reserve_seconds=900,objective_weights='original_base_weights_only',financed_share=.8,experimental_not_adopted=True,no_auto_retry=True,submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
PY
