#!/usr/bin/env bash
#SBATCH --job-name=e5f_evening_gated
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --signal=B:TERM@60
# Six smokes, explicit lead gate, then search/repeats on this SAME allocation.
# Args: STAGE CONTRACT CONTRACT_SHA CONTROLLER CONTROLLER_SHA CUTOFF RUN_ROOT INNER_WRAPPER_SHA
set -euo pipefail
unset NUMBA_DISABLE_JIT APPTAINERENV_NUMBA_DISABLE_JIT SINGULARITYENV_NUMBA_DISABLE_JIT
[ "$#" -eq 8 ] || { echo 'eight explicit pinned arguments required' >&2; exit 2; }
stage=$1; contract=$2; contract_sha=$3; controller=$4; controller_sha=$5
cutoff=$6; run_root=$7; inner_sha=$8
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
[[ "$cutoff" =~ ^[0-9]+$ ]] || exit 2
(( cutoff == 1790558820 && $(date +%s) < cutoff )) || { echo 'search cutoff reached or mismatched' >&2; exit 3; }
[ "${SLURM_CPUS_PER_TASK:-0}" -ge 24 ] || exit 4
case "$run_root" in "$original"/*) physical="$stage/project/${run_root#"$original"/}";; *) echo 'run root outside staged project' >&2; exit 2;; esac
inner="$stage/project/code/cluster/submit_e5f_evening_calibration.sh"
printf '%s  %s\n' "$inner_sha" "$inner" | sha256sum -c -
mkdir "$physical"
child_pid=''
monitor_pid=''
phase() {
  active=$(ps -u "$(id -u)" -o args= | awk -v root="$run_root/" '$1 ~ /\/python/ && index($0,"--stage evaluate") && index($0,root) {n++} END {print n+0}')
  printf '{"phase":"%s","epoch":%s,"job_id":"%s","node":"%s","requested_workers":24,"actual_model_workers":%s}\n' "$1" "$(date +%s)" "$SLURM_JOB_ID" "${SLURMD_NODENAME:-unknown}" "$active" > "$physical/allocation_heartbeat.json.tmp"
  mv "$physical/allocation_heartbeat.json.tmp" "$physical/allocation_heartbeat.json"
}
cleanup() {
  code=$?
  if [ -n "$monitor_pid" ]; then kill -TERM "$monitor_pid" 2>/dev/null || true; fi
  if [ -n "$child_pid" ]; then kill -TERM "$child_pid" 2>/dev/null || true; fi
  if [ "$code" -ne 0 ]; then phase "stopped_exit_$code"; fi
}
trap cleanup EXIT
trap 'exit 143' TERM INT
# Cache reuse never reuses solutions; Numba validates source/CPU compatibility.
if [ ! -e "$stage/numba_cache" ] && [ -d "$stage/native_numba_cache" ]; then
  cp -a "$stage/native_numba_cache" "$stage/numba_cache"
  printf 'native_numba_cache copied remotely; Numba source/CPU validation retained\n' > "$physical/cache_reuse.txt"
fi
run_stage() {
  mode=$1; output=$2; approval=$3; approval_sha=$4
  bash "$inner" "$stage" "$contract" "$contract_sha" "$controller" "$controller_sha" "$cutoff" "$output" "$mode" "$approval" "$approval_sha" &
  child_pid=$!
  (trap - EXIT; trap 'exit 0' TERM INT; while kill -0 "$child_pid" 2>/dev/null; do phase "$mode"; sleep 10; done) &
  monitor_pid=$!
  status=0
  wait "$child_pid" || status=$?
  child_pid=''
  kill -TERM "$monitor_pid" 2>/dev/null || true
  wait "$monitor_pid" 2>/dev/null || true
  monitor_pid=''
  return "$status"
}
phase smoke
run_stage smoke "$run_root/smoke" - -
phase awaiting_explicit_lead_approval
# The lead alone creates these after full smoke/plots AND platform-proof review.
# This wrapper never writes an approval and never infers one from elapsed time.
approval="$physical/lead_search_approval.json"
logical_approval="$run_root/lead_search_approval.json"
sha_file="$physical/lead_search_approval.sha256"
while [ ! -f "$approval" ] || [ ! -f "$sha_file" ]; do
  (( $(date +%s) < cutoff )) || { echo 'no approval before cutoff; stopping' >&2; exit 3; }
  phase awaiting_explicit_lead_approval
  sleep 10
done
(( $(date +%s) < cutoff )) || exit 3
approval_sha=$(tr -d '\r\n' < "$sha_file")
[[ "$approval_sha" =~ ^[0-9a-f]{64}$ ]] || { echo 'invalid approval digest' >&2; exit 2; }
printf '%s  %s\n' "$approval_sha" "$approval" | sha256sum -c -
phase search
run_stage search "$run_root/search" "$logical_approval" "$approval_sha"
phase complete
