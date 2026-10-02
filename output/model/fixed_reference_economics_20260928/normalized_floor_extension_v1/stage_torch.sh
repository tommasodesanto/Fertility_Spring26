#!/usr/bin/env bash
# Stage only; never submits Slurm work.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
remote=/scratch/td2248/projects/normalized_floor_extension_v1
ssh -o BatchMode=yes torch "test -d /scratch/td2248/projects/normalized_calibration_v2/source && test ! -e '$remote' && mkdir -p '$remote/source/output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1' '$remote/logs' '$remote/results'"
scp "$here/run.py" "$here/incumbent.json" "$here/manifest.json" "torch:$remote/source/output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/"
scp "$here/launch_torch.sh" "torch:$remote/launch_torch.sh"
ssh -o BatchMode=yes torch "chmod 755 '$remote/launch_torch.sh' && sha256sum '$remote/source/output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/run.py' '$remote/source/output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/incumbent.json' '$remote/source/output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/manifest.json'"
