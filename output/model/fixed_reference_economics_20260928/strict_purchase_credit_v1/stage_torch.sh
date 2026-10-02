#!/usr/bin/env bash
# Stage the fixed-price diagnostic; do not submit.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
remote=/scratch/td2248/projects/strict_purchase_credit_v1
strict_remote=/scratch/td2248/projects/strict_purchase_sandbox_v1
packet=output/model/fixed_reference_economics_20260928/strict_purchase_credit_v1
ssh -o BatchMode=yes torch "test -d /scratch/td2248/projects/normalized_calibration_v2/source && test -f '$strict_remote/source/code/model/experiments/strict_purchase_sandbox/source/small_credit_lab/engine/household.py' && test ! -e '$remote' && mkdir -p '$remote/source/$packet' '$remote/logs' '$remote/results'"
scp "$here/run.py" "torch:$remote/source/$packet/run.py"
scp -r "$here/input" "torch:$remote/source/$packet/input"
scp "$here/launch_torch.sh" "torch:$remote/launch_torch.sh"
ssh -o BatchMode=yes torch "chmod 755 '$remote/launch_torch.sh' && sha256sum '$remote/source/$packet/run.py' '$remote/source/$packet/input/strict80/collected/run/completed.json' '$strict_remote/source/code/model/experiments/strict_purchase_sandbox/source/small_credit_lab/engine/household.py'"
