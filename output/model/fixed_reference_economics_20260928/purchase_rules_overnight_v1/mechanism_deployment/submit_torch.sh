#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/purchase_mechanism_v1
calibration=/scratch/td2248/projects/purchase_rules_overnight_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate mechanism submission'; exit 2; }
[[ $(date +%s) -lt 1790949600 ]] || { echo '10 AM New York deadline passed'; exit 124; }
/share/apps/anaconda3/2025.06/bin/python - "$remote" "$calibration" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,calibration=map(Path,sys.argv[1:])
manifest=json.loads((remote/'selection/manifest.json').read_text())
assert set(manifest['arms'])=={'hard','quarter'}
for arm,entry in manifest['arms'].items():
 selected=remote/'selection'/f'selected_{arm}.json'
 assert hashlib.sha256(selected.read_bytes()).hexdigest()==entry['selected_json_sha256']
 d=json.loads(selected.read_text())
 assert d['arm']==arm and d['status']=='postchecked' and d['chain']==entry['chain']
 completed=calibration/'results'/f'chain_{entry["chain"]}'/'postcheck/completed.json'
 assert hashlib.sha256(completed.read_bytes()).hexdigest()==entry['completed_sha256']
 assert json.loads(completed.read_text())['status']=='selected_numerically_verified'
 assert (completed.parent/'selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz').is_file()
PY
smoke=$(sbatch --parsable --array=0-1%2 --time=00:30:00 --export=ALL,MECHANISM_SMOKE=1 launch_torch.sh)
smoke=${smoke%%;*}
production=$(sbatch --parsable --dependency="afterok:$smoke" --array=0-11%12 --time=04:00:00 --export=ALL,MECHANISM_SMOKE=0 launch_torch.sh)
production=${production%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$smoke" "$production" <<'PY'
import json,sys,time
from pathlib import Path
r=dict(status='submitted_native_smoke_dependency_gated',smoke_array_job_id=sys.argv[1],production_array_job_id=sys.argv[2],smoke_cases=2,production_cases=12,dependency='afterok:'+sys.argv[1],cpus_per_case=1,memory_GiB_per_case=32,threads_per_case=1,wall_seconds_per_case=14400,maximum_native_policy_calls_per_case=1024,absolute_deadline_epoch=1790949600,no_auto_retry=True,submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
PY
