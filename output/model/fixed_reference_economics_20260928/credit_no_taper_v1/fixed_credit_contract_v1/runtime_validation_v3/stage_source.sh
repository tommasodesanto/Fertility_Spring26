#!/usr/bin/env bash
set -euo pipefail
local_source="$(cd "$(dirname "$0")" && pwd)"
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
overlay="$frozen/output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/overlay"
remote=/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v3
[[ -f "$local_source/run_runtime_validation.py" && -f "$overlay/solver.py" ]] || { echo "uploaded packet or frozen overlay absent" >&2; exit 2; }
[[ ! -e "$remote/source" ]] || { echo "refusing to overwrite staged source" >&2; exit 2; }
mkdir -p "$remote/source/overlay" "$remote/results" "$remote/cache"
cp "$local_source/run_runtime_validation.py" "$local_source/plan.json" "$local_source/launch_runtime_validation.sh" "$local_source/launch_smoke.sh" "$local_source/launch_control.sh" "$remote/source/"
cp "$overlay/parameters.py" "$overlay/solver.py" "$overlay/kernels.py" "$remote/source/overlay/"
chmod +x "$remote/source/launch_runtime_validation.sh" "$remote/source/launch_smoke.sh" "$remote/source/launch_control.sh"
sha256sum "$remote/source/run_runtime_validation.py" "$remote/source/plan.json" "$remote/source/launch_runtime_validation.sh" "$remote/source/launch_smoke.sh" "$remote/source/launch_control.sh" "$remote/source/overlay/"*.py > "$remote/source/source.sha256"
chmod -R a-w "$remote/source"
echo "staged read-only source at $remote/source"
