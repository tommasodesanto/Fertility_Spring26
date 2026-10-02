#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a
selection_store=/scratch/td2248/projects/purchase_mechanism_v1
cd "$remote"
[[ ! -e submission_receipt.json ]] || { echo 'Refusing duplicate mechanism submission'; exit 2; }
[[ $(date +%s) -lt 1790949600 ]] || { echo '10 AM New York deadline passed'; exit 124; }
/share/apps/anaconda3/2025.06/bin/python - "$remote" "$selection_store" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote=Path(sys.argv[1]);selection_store=Path(sys.argv[2])
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''): h.update(block)
 return h.hexdigest()
assert hashlib.sha256((remote/'inventory.json').read_bytes()).hexdigest()=='7cd28f68e96a18b24b3bea9e9b709038576444de20032a9538a82c37574247a3'
route=json.loads((remote/'reviewed_routing_receipt.json').read_text())
assert route['launch_sha256']==sha(remote/'launch_torch_reviewed.sh')
assert route['submit_sha256']==sha(remote/'submit_torch_reviewed.sh')
manifest_path=selection_store/'selection/manifest.json'
assert hashlib.sha256(manifest_path.read_bytes()).hexdigest()=='19db8dfcc928c4cdea5810f70c4e62a9de65372bd038646d0bc2d49532913879'
manifest=json.loads(manifest_path.read_text())
assert set(manifest['arms'])=={'hard','quarter'}
for arm,entry in manifest['arms'].items():
 selected=selection_store/'selection'/f'selected_{arm}.json'
 assert hashlib.sha256(selected.read_bytes()).hexdigest()==entry['selected_json_sha256']
 d=json.loads(selected.read_text())
 assert d['arm']==arm and d['status']=='postchecked' and d['chain']==entry['chain']
 snapshot=selection_store/'selected_postchecks'/f'chain_{entry["chain"]}'
 assert entry['snapshot_remote_root']==str(snapshot) and d['snapshot_remote_root']==str(snapshot)
 assert d['origin']==entry['origin'] and d['remote_root']==entry['physical_remote_root']
 assert d.get('parent_remote_root')==entry['parent_remote_root']
 completed=snapshot/'postcheck/completed.json'
 assert hashlib.sha256(completed.read_bytes()).hexdigest()==entry['completed_sha256']
 assert json.loads(completed.read_text())['status']=='selected_numerically_verified'
 report=completed.parent/'selected_postcheck/phase_b_ge/selected_root'
 for name,digest in entry['report_sha256'].items():
  assert hashlib.sha256((report/name).read_bytes()).hexdigest()==digest
 arrays=completed.parent/'selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz'
 assert sha(arrays)==entry['native_arrays_sha256']
PY
smoke=$(sbatch --parsable --array=0-1%2 --time=00:30:00 --export=ALL,MECHANISM_SMOKE=1 launch_torch_reviewed.sh)
smoke=${smoke%%;*}
production=$(sbatch --parsable --dependency="afterok:$smoke" --array=0-11%12 --time=04:00:00 --export=ALL,MECHANISM_SMOKE=0 launch_torch_reviewed.sh)
production=${production%%;*}
/share/apps/anaconda3/2025.06/bin/python - "$smoke" "$production" <<'PY'
import json,sys,time
from pathlib import Path
r=dict(source_root='/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a',selection_store='/scratch/td2248/projects/purchase_mechanism_v1',status='submitted_native_smoke_dependency_gated',smoke_array_job_id=sys.argv[1],production_array_job_id=sys.argv[2],smoke_cases=2,production_cases=12,dependency='afterok:'+sys.argv[1],cpus_per_case=1,memory_GiB_per_case=32,threads_per_case=1,wall_seconds_per_case=14400,maximum_native_policy_calls_per_case=1024,absolute_deadline_epoch=1790949600,no_auto_retry=True,submission_epoch=time.time())
Path('submission_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
PY
