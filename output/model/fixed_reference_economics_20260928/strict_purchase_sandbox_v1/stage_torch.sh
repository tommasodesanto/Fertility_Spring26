#!/usr/bin/env bash
# Stage the isolated source and one-case runner; never submit.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
remote=/scratch/td2248/projects/strict_purchase_sandbox_v1
sandbox="$repo/code/model/experiments/strict_purchase_sandbox"
ssh -o BatchMode=yes torch "test -d /scratch/td2248/projects/normalized_calibration_v2/source && test ! -e '$remote' && mkdir -p '$remote/source/output/model/fixed_reference_economics_20260928/strict_purchase_sandbox_v1' '$remote/source/code/model/experiments' '$remote/logs' '$remote/results'"
scp "$here/run.py" "$here/incumbent.json" "$here/manifest.json" "torch:$remote/source/output/model/fixed_reference_economics_20260928/strict_purchase_sandbox_v1/"
scp -r "$sandbox" "torch:$remote/source/code/model/experiments/"
scp "$here/launch_torch.sh" "torch:$remote/launch_torch.sh"
ssh -o BatchMode=yes torch "chmod 755 '$remote/launch_torch.sh' && sha256sum '$remote/source/output/model/fixed_reference_economics_20260928/strict_purchase_sandbox_v1/run.py' '$remote/source/output/model/fixed_reference_economics_20260928/strict_purchase_sandbox_v1/manifest.json' '$remote/source/code/model/experiments/strict_purchase_sandbox/source/small_credit_lab/engine/household.py' '$remote/source/code/model/experiments/strict_purchase_sandbox/source/small_credit_lab/engine/kernels.py'"
