#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
remote=/scratch/td2248/projects/purchase_fresh_calibration_v1
"${PYTHON:-python3}" - "$here" <<'PY'
import hashlib,json,sys
from pathlib import Path
here=Path(sys.argv[1]);names=('design.json','driver.py','launch_torch.sh','preflight_torch.sh','test_design.py')
pins={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in names}
(here/'source_sha256.json').write_text(json.dumps(pins,sort_keys=True,indent=2)+'\n')
print(json.dumps(pins,sort_keys=True))
PY
ssh torch "mkdir -p '$remote/source' '$remote/results' '$remote/logs' '$remote/preflight'"
scp "$here/design.json" "$here/driver.py" "$here/launch_torch.sh" "$here/preflight_torch.sh" "$here/test_design.py" torch:"$remote/source/"
scp "$here/launch_torch.sh" "$here/preflight_torch.sh" "$here/source_sha256.json" torch:"$remote/"
ssh torch "'/share/apps/anaconda3/2025.06/bin/python' - '$remote'" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);pins=json.loads((root/'source_sha256.json').read_text())
for name,digest in pins.items():
 p=root/'source'/name
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
for name in ('launch_torch.sh','preflight_torch.sh'):
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==pins[name]
print(json.dumps(dict(status='staged_source_only',source_sha256=pins)))
PY
