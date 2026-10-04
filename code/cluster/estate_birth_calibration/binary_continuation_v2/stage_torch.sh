#!/usr/bin/env bash
# Build and stage sources, then derive authenticated starts. Does not submit jobs.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
remote=/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2
recovery=/scratch/td2248/projects/estate_birth_recovery_20261004_v1
parent=/scratch/td2248/projects/estate_birth_calibration_20261003_v3
python3 "$here/build_stage.py"
local="$repo/output/model/experiments/birth_count_choice/estate_a_binary_continuation_20261004_v2/deployment"
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir -p '$remote/logs' '$remote/control'"
scp "$local/stage.tar.gz" torch:"$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && chmod 755 ./*.sh && /share/apps/anaconda3/2025.06/bin/python verify_stage.py --host && /share/apps/anaconda3/2025.06/bin/python prepare_starts.py --recovery '$recovery' --parent '$parent' --out '$remote/control' && /share/apps/anaconda3/2025.06/bin/python verify_starts.py '$remote/control/starts.json'"
scp "torch:$remote/control/starts.json" "$local/starts.json"
scp "torch:$remote/control/starts_receipt.json" "$local/starts_receipt.json"
python3 "$here/verify_starts.py" "$local/starts.json" --structure-only
