#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
remote=/scratch/td2248/projects/purchase_restart_controller_v2
"${PYTHON:-python3}" - "$here" <<'PY'
import hashlib,json,sys
from pathlib import Path
here=Path(sys.argv[1]);names=('controller.py','test_controller.py','launch_torch.sh')
pins={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in names}
(here/'source_sha256.json').write_text(json.dumps(pins,sort_keys=True,indent=2)+'\n')
print(json.dumps(pins,sort_keys=True))
PY
ssh torch "mkdir -p '$remote/source' '$remote/results' '$remote/logs'"
scp "$here/controller.py" "$here/test_controller.py" "$here/launch_torch.sh" torch:"$remote/source/"
scp "$here/launch_torch.sh" "$here/source_sha256.json" torch:"$remote/"
ssh torch "'/share/apps/anaconda3/2025.06/bin/python' - '$remote'" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);pins=json.loads((root/'source_sha256.json').read_text())
for name,digest in pins.items():
 p=root/'source'/name
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
assert hashlib.sha256((root/'launch_torch.sh').read_bytes()).hexdigest()==pins['launch_torch.sh']
print(json.dumps(dict(status='staged_source_only',source_sha256=pins)))
PY
