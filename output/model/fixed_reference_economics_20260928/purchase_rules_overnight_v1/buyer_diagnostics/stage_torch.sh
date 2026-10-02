#!/usr/bin/env bash
# Stage only the buyer addendum; never alter calibration or mechanism pins.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
remote=/scratch/td2248/projects/purchase_buyer_diagnostics_v4
python="$here/../../../../../code/model/.venv/bin/python"
"$python" "$here/build_stage.py" > "$here/stage_build_receipt.json"
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir '$remote'"
scp "$here/buyer_diagnostics_stage_v4.tar.gz" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && mkdir -p logs results && chmod 755 source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/buyer_diagnostics/launch_readout.sh && sha256sum stage.tar.gz"
