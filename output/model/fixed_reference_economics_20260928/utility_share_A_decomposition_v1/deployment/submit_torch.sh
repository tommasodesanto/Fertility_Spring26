#!/usr/bin/env bash
set -euo pipefail
cd /scratch/td2248/projects/utility_share_A_decomposition_v1
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate submission'; exit 2; }
[[ -f GO_REVIEWED ]] || { echo 'Requires lead GO_REVIEWED'; exit 2; }
job=$(sbatch --parsable launch_torch.sh)
job=${job%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$job" <<'PY'
import json,sys,time
from pathlib import Path
r=dict(job_id=sys.argv[1],submitted_epoch=time.time(),cpus=1,memory_GiB=24,threads=1,total_wall_seconds=1200,lifecycle_cap=3,per_cell_seconds=600,no_auto_retry=True,experimental_not_adopted=True)
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
PY
