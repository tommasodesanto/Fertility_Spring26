#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_winner31_credit_ge_replacement_v1
fixed=/scratch/td2248/projects/utility_floor_winner31_responses_v1
q0="$fixed/results/responses/04_lifetime_repayment_only_p1.00"
cd "$remote"
[[ -f "$q0/receipt.json" && -f "$q0/closure.json" ]] || { echo 'Exact winner q0 expanded-credit result missing'; exit 2; }
(( $(date +%s) < 1790870869 )) || { echo 'Original absolute deadline reached'; exit 124; }
[[ ! -e submission_receipt.json && ! -e results ]] || { echo 'Refusing duplicate submission'; exit 2; }
job=$(sbatch --parsable --time=00:40:00 launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY_SUBMIT'
import json,sys,time
from pathlib import Path
r=dict(status='submitted_replacement_after_presolve_driver_failure',job_id=sys.argv[1],submission_epoch=time.time(),original_absolute_deadline_epoch=1790870869,original_failed_job_id='18958106',prior_attempts_reserved=1,cpu=1,memory_gib=24,threads=1,maximum_new_lifecycle_solves=11,case_seconds=300,total_seconds_from_original_start=2400,no_auto_retry=True,natural_support_certified=False)
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
PY_SUBMIT
