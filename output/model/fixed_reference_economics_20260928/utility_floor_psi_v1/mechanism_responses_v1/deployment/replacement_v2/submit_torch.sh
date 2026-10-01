#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_mechanism_responses_v1_replacement2
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate submission'; exit 2; }
submission_epoch=$(date +%s)
deadline_epoch=1790833605
remaining=$((deadline_epoch-submission_epoch))
[[ "$remaining" -gt 600 ]] || { echo "Insufficient original deadline remaining"; exit 2; }
minutes=$(((remaining+59)/60))
printf '%s\n' "$deadline_epoch" > common_deadline_epoch
job=$(sbatch --parsable --time="$minutes" launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" "$submission_epoch" "$deadline_epoch" <<'PY'
import json,sys
from pathlib import Path
r=dict(status='submitted',job_id=sys.argv[1],submission_epoch=int(sys.argv[2]),deadline_epoch=int(sys.argv[3]),wall_seconds=int(sys.argv[3])-int(sys.argv[2]),original_failed_job='18922355',prior_completed_lifecycle_job='18925514',prior_completed_lifecycle_solves=1,cumulative_mechanism_maximum_lifecycle_solves=7,original_deadline_preserved=True,cpu=1,memory_gib=16,threads=1,maximum_lifecycle_solves=6,case_seconds=600,no_auto_retry=True,natural_support_certified=False,credit_ge_launched=False)
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
PY
