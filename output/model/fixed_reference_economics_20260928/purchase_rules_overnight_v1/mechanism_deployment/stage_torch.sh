#!/usr/bin/env bash
# Source-only staging. Does not create or submit model jobs.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
remote=/scratch/td2248/projects/purchase_mechanism_v1
python3 "$here/build_stage.py" --archive "$here/mechanism_stage_v2.tar.gz"
archive="$here/mechanism_stage_v2.tar.gz"
ssh -o BatchMode=yes torch "test ! -e '$remote/submission_receipt.json' && mkdir -p '$remote/logs' '$remote/results' '$remote/selection'"
scp "$archive" "torch:$remote/stage_v2.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage_v2.tar.gz && chmod 755 launch_torch.sh submit_torch.sh && sha256sum stage_v2.tar.gz"
