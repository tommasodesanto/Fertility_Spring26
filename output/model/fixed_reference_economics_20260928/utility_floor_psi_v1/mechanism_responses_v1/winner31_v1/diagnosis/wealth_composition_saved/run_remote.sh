#!/usr/bin/env bash
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_winner31_responses_v1
analysis="$remote/analysis_winner31_wealth"
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
test -f "$remote/results/responses/q0_reference_inherited_states.npz"
test -f "$remote/results/responses/00_reference_p1.00/solution_arrays.npz"
test -f "$remote/results/responses/04_lifetime_repayment_only_p1.00/solution_arrays.npz"
test ! -e "$analysis/result"
mkdir -p "$analysis/numba_cache" "$analysis/matplotlib"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR=/work/analysis/numba_cache MPLCONFIGDIR=/work/analysis/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$remote/results:/work/results:ro" --bind "$analysis:/work/analysis:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$remote/source/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1:$repo/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1:ro")
while IFS= read -r rel; do binds+=(--bind "$remote/source/$rel:$repo/$rel:ro"); done < "$remote/mounts.txt"
start=$(date +%s)
set +e
timeout --signal=TERM --kill-after=5s 180s apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/analysis/reduce.py --packet "$repo/$packet" --root /work/results/responses --out /work/analysis/result > "$analysis/stdout.txt" 2> "$analysis/stderr.txt"
code=$?
set -e
printf 'exit=%s elapsed=%s\n' "$code" "$(( $(date +%s)-start ))"
exit "$code"
