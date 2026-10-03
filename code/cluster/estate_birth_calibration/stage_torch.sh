#!/usr/bin/env bash
# Stage new immutable snapshot and initialize without solves. Never submit here.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
remote=/scratch/td2248/projects/estate_birth_calibration_20261003_v2
python3 "$here/build_stage.py"
archive="$repo/output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt2/stage.tar.gz"
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir -p '$remote/logs'"
scp "$archive" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && chmod 755 ./*.sh && /share/apps/anaconda3/2025.06/bin/python verify_stage.py --host && ./launch_torch.sh preflight && PREFLIGHT_TASK=5 ./launch_torch.sh preflight"
