#!/usr/bin/env bash
set -euo pipefail
: "${FLOOR_GATE_JOB:?Provide the successful external gate job ID}"
[[ "$FLOOR_GATE_JOB" =~ ^[0-9]+$ ]] || exit 2
remote=/scratch/td2248/projects/normalized_floor_calibration_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate submission'; exit 2; }
# The launcher authenticates complete gate evidence again before any solve.
state=$(sacct -n -X -j "$FLOOR_GATE_JOB" --format=State --parsable2 | head -1)
[[ "$state" == COMPLETED ]] || { echo 'External gate must already be completed successfully'; exit 2; }
prod=$(sbatch --parsable --dependency="afterok:$FLOOR_GATE_JOB" --export="ALL,FLOOR_GATE_JOB=$FLOOR_GATE_JOB" --array=0-23%24 --time=03:00:00 launch_torch.sh)
prod=${prod%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$FLOOR_GATE_JOB" "$prod" <<'PY'
import json,sys,time
from pathlib import Path
r=dict(status='submitted_external_dependency_gated',gate_job_id=sys.argv[1],array_job_id=sys.argv[2],chains=24,dependency=f'afterok:{sys.argv[1]}',cpus_per_chain=1,memory_GiB_per_chain=24,threads_per_chain=1,actual_start_seconds_per_chain=10800,maximum_objective_calls_per_chain=150,final_native_reserve_seconds=900,objective_weights='original_base_weights_only',experimental_not_adopted=True,no_auto_retry=True,submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
PY
