#!/usr/bin/env bash
# Run through ordinary authorized Torch SSH. Preparation only; no Python/model calls.
set -euo pipefail
old=/scratch/td2248/projects/fertility_evening_calibration_20260927_v1
new=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
[ -d "$old/project" ]
[ ! -e "$new" ] || { echo 'Fresh stage required; inspect existing attempt' >&2; exit 2; }
mkdir "$new"
# Physical copy: old results stay immutable; omit the two large case trees.
rsync -a --exclude='/output/model/evening_calibration_20260927/gated_v1/' \
  --exclude='/output/model/evening_calibration_20260927/gated_v2/' \
  --exclude='__pycache__/' "$old/project/" "$new/project/"
rel=output/model/evening_calibration_20260927/gated_v2/search
mkdir -p "$new/project/$rel/selected_export"
for name in initial_0347_block repeat_0364_block repeat_0365_block; do
  cp -a "$old/project/$rel/$name" "$new/project/$rel/$name"
done
# Materialize the selected export's checkpoint symlink against the new copied case.
cp -a "$old/project/$rel/selected_export/block" "$new/project/$rel/selected_export/block"
# Exact-root absolute symlink remains valid inside the explicit Apptainer bind.
# Cache is non-scientific and Numba validates source/CPU compatibility.
if [ -d "$old/numba_cache" ]; then cp -a "$old/numba_cache" "$new/numba_cache"; fi
printf 'old=%s\nnew=%s\nstatus=remote_copy_complete_no_imports_no_solves\n' "$old" "$new" > "$new/staging_receipt.txt"
