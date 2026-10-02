#!/usr/bin/env bash
# Stage two-case quarter-saving diagnostic; do not submit.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
remote=/scratch/td2248/projects/quarter_saving_benchmarks_v1
packet=output/model/fixed_reference_economics_20260928/quarter_saving_benchmarks_v1
sandbox="$repo/code/model/experiments/quarter_saving_sandbox"
ssh -o BatchMode=yes torch "test -d /scratch/td2248/projects/normalized_calibration_v2/source && test ! -e '$remote' && mkdir -p '$remote/source/$packet' '$remote/source/code/model/experiments' '$remote/logs' '$remote/results'"
scp "$here/run.py" "$here/manifest.json" "torch:$remote/source/$packet/"
scp -r "$here/input" "torch:$remote/source/$packet/input"
scp -r "$sandbox" "torch:$remote/source/code/model/experiments/"
scp "$here/launch_torch.sh" "torch:$remote/launch_torch.sh"
ssh -o BatchMode=yes torch "chmod 755 '$remote/launch_torch.sh' && sha256sum '$remote/source/$packet/run.py' '$remote/source/$packet/manifest.json' '$remote/source/code/model/experiments/quarter_saving_sandbox/source/small_credit_lab/engine/household.py' '$remote/source/code/model/experiments/quarter_saving_sandbox/source/small_credit_lab/engine/kernels.py'"
