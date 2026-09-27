#!/usr/bin/env bash
#SBATCH --job-name=e5f_evening
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cpu_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --signal=B:TERM@60
# Lead submits only after the exact-loop smoke/source/plan review.
# Usage: sbatch THIS STAGE CONTRACT CONTRACT_SHA CONTROLLER CONTROLLER_SHA SEARCH_CUTOFF OUTPUT MODE APPROVAL APPROVAL_SHA
set -euo pipefail
if [ "$#" -ne 10 ]; then
  echo 'usage: STAGE CONTRACT CONTRACT_SHA CONTROLLER CONTROLLER_SHA SEARCH_CUTOFF OUTPUT MODE APPROVAL APPROVAL_SHA' >&2
  exit 2
fi
stage=$1
plan=$2
plan_sha=$3
controller=$4
controller_sha=$5
search_cutoff=$6
output=$7
mode=$8
approval=$9
approval_sha=${10}
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
hard_end=1790561520 # 2026-09-28 02:12 UTC; no queue-time extension.
now=$(date +%s)
[[ "$search_cutoff" =~ ^[0-9]+$ ]] || { echo 'integer search cutoff required' >&2; exit 2; }
if (( now >= search_cutoff || now >= hard_end || search_cutoff > hard_end )); then
  echo 'expired search window; refusing late numerical start' >&2
  exit 3
fi
if [ "${SLURM_CPUS_PER_TASK:-0}" -lt 24 ]; then
  echo '24 allocated CPUs required; refuse unsupported resource allocation' >&2
  exit 4
fi
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$stage/numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"
# Resolve original absolute paths only through the explicit staged-root mapping.
case "$plan" in "$original"/*) staged_plan="$stage/project/${plan#"$original"/}";; *) echo 'plan outside authenticated project' >&2; exit 2;; esac
case "$controller" in "$original"/*) staged_controller="$stage/project/${controller#"$original"/}";; *) echo 'controller outside authenticated project' >&2; exit 2;; esac
printf '%s  %s\n' "$plan_sha" "$staged_plan" | sha256sum -c -
printf '%s  %s\n' "$controller_sha" "$staged_controller" | sha256sum -c -
args=(--stage "$mode" --contract "$plan" --output "$output")
case "$mode" in
  prepare|smoke) ;;
  search)
    case "$approval" in "$original"/*) staged_approval="$stage/project/${approval#"$original"/}";; *) staged_approval="$approval";; esac
    printf '%s  %s\n' "$approval_sha" "$staged_approval" | sha256sum -c -
    args+=(--approval "$approval" --approval-sha256 "$approval_sha")
    ;;
  *) echo 'mode must be prepare, smoke or search' >&2; exit 2;;
esac
export EXPECTED_E5F_EVENING_SHA256="$plan_sha"
remaining=$((hard_end - $(date +%s) - 10)) # reserve hard-kill grace inside end
(( remaining > 0 )) || exit 3
exec timeout --signal=TERM --kill-after=10 "$remaining" \
  apptainer exec --bind "$stage/project:$original" --pwd "$original" \
  /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python \
  "$controller" "${args[@]}"
