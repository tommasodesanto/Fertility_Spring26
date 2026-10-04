#!/usr/bin/env bash
# Explicit lead action: build, authenticate, and stage; never submits Slurm.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); repo=$(cd "$here/../../.." && pwd)
remote=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3
python="$repo/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python"
deployment="$repo/output/model/experiments/ces_normalized_shares/overnight_v1/deployment/attempt3"
archive="$deployment/stage.tar.gz"
receipt="$deployment/stage_receipt.json"
# Reuse authenticated bytes after an interrupted transfer; never rebuild a
# retained attempt merely because the network disconnected.
if [[ ! -f "$archive" || ! -f "$receipt" ]]; then
  "$python" "$here/build_stage.py"
fi
expected_sha=$("$python" - "$archive" "$receipt" <<'PYTHON'
import hashlib, json, sys
from pathlib import Path
archive, receipt = map(Path, sys.argv[1:])
expected = json.loads(receipt.read_text())["sha256"]
actual = hashlib.sha256(archive.read_bytes()).hexdigest()
if actual != expected:
    raise SystemExit("Retained stage archive differs from its receipt; refusing transfer")
print(actual)
PYTHON
)
ssh -o BatchMode=yes torch "test ! -e '$remote/inventory.json' && mkdir -p '$remote/logs'"
rsync --partial --append -e 'ssh -o BatchMode=yes' "$archive" "torch:$remote/stage.tar.gz"
# rsync --append is resumable on the Mac, but authentication is a separate
# full remote digest check before any archive bytes are extracted or executed.
ssh -o BatchMode=yes torch "bash -s -- '$remote' '$expected_sha'" <<'REMOTE'
set -euo pipefail
remote=$1
expected_sha=$2
cd "$remote"
test ! -e inventory.json
actual_sha=$(sha256sum stage.tar.gz)
actual_sha=${actual_sha%% *}
if [[ "$actual_sha" != "$expected_sha" ]]; then
  echo "Remote stage archive differs from its receipt; refusing extraction" >&2
  exit 1
fi
tar -xzf stage.tar.gz
chmod 755 ./*.sh
/share/apps/anaconda3/2025.06/bin/python verify_stage.py --host
REMOTE
ssh -o BatchMode=yes torch "cd '$remote' && ./launch_torch.sh preflight"
