#!/usr/bin/env bash
# Replace only this isolated stage after the corrected archive upload is approved.
# The failed first archive remains at stage_attempt1.tar.gz.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
remote=/scratch/td2248/projects/soft_timing_calibration_20261002_v1
python3 "$here/build_stage.py"
archive="$repo/output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/stage.tar.gz"
ssh -o BatchMode=yes torch "test -d '$remote' && test -e '$remote/stage_attempt1.tar.gz' && test ! -e '$remote/submission_receipt.json' && test ! -e '$remote/results' && test ! -e '$remote/stage.tar.gz'"
scp "$archive" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && chmod 755 ./*.sh && /share/apps/anaconda3/2025.06/bin/python verify_stage.py --host && sha256sum stage.tar.gz && ./launch_torch.sh preflight"
