#!/usr/bin/env bash
# Source-only package staging. This script never submits jobs.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../../.." && pwd)
remote=/scratch/td2248/projects/purchase_rules_overnight_v1
python3 "$here/build_stage.py"
archive="$here/purchase_rules_overnight_v1_stage.tar.gz"
ssh -o BatchMode=yes torch "mkdir -p '$remote/logs' '$remote/results' && test ! -e '$remote/submission_receipt.json'"
scp "$archive" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && chmod 755 launch_torch.sh submit_torch.sh && sha256sum stage.tar.gz"
