#!/usr/bin/env bash
# Torch-only Slurm payload.  This script never calls sbatch and never retries.
# The lead stages the controller/source and invokes sbatch explicitly after review.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
set -euo pipefail

stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
overlay=/scratch/td2248/projects/fixed_reference_transition_20260928
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
relative=output/model/fixed_reference_transition_20260928
controller_rel="${E5F_PE_CONTROLLER_RELATIVE:?set lead-staged controller relative path}"
source_rel="${E5F_PE_SOURCE_RELATIVE:?set lead-staged estimator source relative path}"
plan_rel="${E5F_PE_PLAN_RELATIVE:?set enabled plan relative to project root}"
output_rel="${E5F_PE_OUTPUT_RELATIVE:?set new output relative to project root}"
plan_sha="${E5F_PE_PLAN_SHA256:?set pinned plan SHA-256}"
controller_sha="${E5F_PE_CONTROLLER_SHA256:?set controller SHA-256}"
launcher_sha="${E5F_PE_LAUNCHER_SHA256:?set launcher SHA-256}"
expected_id="${EXPECTED_E5F_IDENTIFICATION_SHA256:-}"

require_overlay_relative() {
  local value="$1" physical
  case "$value" in
    /*|*'//'*) return 1 ;;
  esac
  IFS=/ read -r -a parts <<< "$value"
  for part in "${parts[@]}"; do
    [[ -n "$part" && "$part" != . && "$part" != .. ]] || return 1
  done
  [[ "$value" == "$relative/"* ]] || return 1
  physical="$overlay/${value#"$relative/"}"
  printf '%s\n' "$physical"
}

test "$expected_id" = 68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
controller_physical="$(require_overlay_relative "$controller_rel")" || { echo 'controller path is outside overlay' >&2; exit 2; }
source_physical="$(require_overlay_relative "$source_rel")" || { echo 'source path is outside overlay' >&2; exit 2; }
plan_physical="$(require_overlay_relative "$plan_rel")" || { echo 'plan path is outside overlay' >&2; exit 2; }
output_physical="$(require_overlay_relative "$output_rel")" || { echo 'output path is outside overlay' >&2; exit 2; }
test -r "$controller_physical"
test -d "$source_physical"
test -r "$plan_physical"
test "$(sha256sum "$controller_physical" | awk '{print $1}')" = "$controller_sha"
test "$(sha256sum "$0" | awk '{print $1}')" = "$launcher_sha"
test "$(sha256sum "$plan_physical" | awk '{print $1}')" = "$plan_sha"
test ! -e "$output_physical"

module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export NUMBA_CACHE_DIR="$overlay/four_shock_v1/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
exec apptainer exec --bind "$stage/project:$original:ro" --bind "$overlay:$original/$relative:rw" --pwd "$original" \
  /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$controller_rel" execute \
  --python /share/apps/anaconda3/2025.06/bin/python --plan "$plan_rel" --plan-sha256 "$plan_sha" \
  --source-dir "$source_rel" --output "$output_rel"
