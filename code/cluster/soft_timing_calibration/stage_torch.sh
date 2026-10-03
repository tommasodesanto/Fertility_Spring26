#!/usr/bin/env bash
# Stage immutable source and run zero-solve verification. Never submits Slurm.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
remote=/scratch/td2248/projects/soft_timing_calibration_20261002_v1
python3 "$here/build_stage.py"
archive="$repo/output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/stage.tar.gz"
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir -p '$remote'"
scp "$archive" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && chmod 755 ./*.sh && /share/apps/anaconda3/2025.06/bin/python verify_stage.py --host && sha256sum stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && ./launch_torch.sh preflight"
