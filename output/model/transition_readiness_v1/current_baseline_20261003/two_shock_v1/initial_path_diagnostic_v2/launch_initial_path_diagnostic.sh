#!/usr/bin/env bash
#SBATCH --job-name=estate_path_diagnostic
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cl
#SBATCH --output=/scratch/td2248/projects/current_estate_initial_path_20261004_v2/logs/%x-%j.out
set -euo pipefail

repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
local_addon="$repo/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/initial_path_diagnostic_v2"
remote_addon=/scratch/td2248/projects/current_estate_initial_path_20261004_v2
local_v5="$repo/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5"
remote_v5=/scratch/td2248/projects/current_estate_two_shock_20261004_v5
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
expected_config_sha=${1:-}

[[ "$#" == 1 ]] || { echo 'Exactly one config SHA256 argument is required'; exit 2; }
[[ "$expected_config_sha" =~ ^[[:xdigit:]]{64}$ ]] || { echo 'Pass the exact 64-character config SHA256 as the only argument'; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'A single Slurm job is required'; exit 2; }
[[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'Arrays forbidden'; exit 2; }
[[ -z "${SLURM_NTASKS:-}" || "$SLURM_NTASKS" == 1 ]] || { echo 'Exactly one Slurm task required'; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-1}" == 1 ]] || { echo 'Exactly one CPU required'; exit 2; }
[[ "${SLURM_MEM_PER_NODE:-24576}" -le 24576 ]] || { echo '24GiB memory cap exceeded'; exit 2; }

module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH

mkdir -p "$remote_addon/jobs" "$remote_addon/results" "$remote_addon/logs"
claim="$remote_addon/jobs/mapping.claim"
mkdir "$claim" || { echo "Refusing duplicate or unknown mapping outcome: $claim"; exit 2; }
printf '%s\n' "${SLURM_JOB_ID}" > "$claim/slurm_job_id"
printf '%s\n' "$expected_config_sha" > "$claim/config_sha256"
job="$remote_addon/jobs/mapping_${SLURM_JOB_ID}"
out="$remote_addon/results/map_v1"
mkdir "$job" || { echo "Refusing existing or unknown job: $job"; exit 2; }
mkdir "$out" || { echo "Refusing existing or unknown mapping output: $out"; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
visible_out="$local_addon/results/map_v1"

start_epoch=$(date +%s); deadline_epoch=$((start_epoch+1800))
child_pid=''; heartbeat_pid=''
stop_child() {
 if [[ -n "$child_pid" ]] && kill -0 "$child_pid" 2>/dev/null; then
  kill -TERM -- "-$child_pid" 2>/dev/null || true
  for _ in {1..10}; do kill -0 "$child_pid" 2>/dev/null || break; sleep 1; done
  kill -KILL -- "-$child_pid" 2>/dev/null || true
  wait "$child_pid" 2>/dev/null || true
 fi
}
terminal() {
 status=$?; trap - EXIT TERM INT
 stop_child
 [[ -z "$heartbeat_pid" ]] || kill "$heartbeat_pid" 2>/dev/null || true
 "$python" - "$job/launcher_terminal.json" "$status" "$start_epoch" "$deadline_epoch" <<'PYRECEIPT'
import json,os,sys,time
from pathlib import Path
path,status,start,deadline=sys.argv[1:]
Path(path).write_text(json.dumps(dict(status='complete' if int(status)==0 else 'failed',exit_code=int(status),start_epoch=int(start),deadline_epoch=int(deadline),finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),sort_keys=True,indent=2)+'\n')
PYRECEIPT
 cp "$job/launcher_terminal.json" "$out/launcher_terminal.json"
 exit "$status"
}
trap terminal EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

config="$remote_addon/config.json"
actual_config_sha=$(sha256sum "$config" | awk '{print $1}')
[[ "$actual_config_sha" == "$expected_config_sha" ]] || { echo 'Remote config SHA256 differs from the supplied pin'; exit 2; }
[[ -f "$remote_addon/diagnose_initial_path.py" ]] || { echo 'Diagnostic driver is missing'; exit 2; }

"$python" - "$job/launcher_start.json" "$start_epoch" "$deadline_epoch" "$actual_config_sha" <<'PYRECEIPT'
import json,os,sys
from pathlib import Path
path,start,deadline,sha=sys.argv[1:]
Path(path).write_text(json.dumps(dict(start_epoch=int(start),deadline_epoch=int(deadline),config_sha256=sha,output='results/map_v1/run',cpus=1,numba_threads=1,blas_threads=1,memory_gib=24,slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),sort_keys=True,indent=2)+'\n')
PYRECEIPT
cp "$job/launcher_start.json" "$out/launcher_start.json"
(while true; do date -u +%Y-%m-%dT%H:%M:%SZ > "$job/launcher_heartbeat.txt"; sleep 60; done) >/dev/null 2>&1 < /dev/null &
heartbeat_pid=$!
run_bounded() {
 local remaining=$((deadline_epoch-$(date +%s)-15))
 [[ "$remaining" -gt 0 ]] || return 124
 setsid timeout --signal=TERM --kill-after=10s "${remaining}s" "$@" <&0 & child_pid=$!
 wait "$child_pid"; local status=$?; child_pid=''; return "$status"
}

# Verify the authenticated immutable runtime on the host and in the container before the diagnostic.
run_bounded "$python" "$remote_v5/prepare_two_shock_torch.py" verify --stage "$remote_v5" > "$out/host_verification.json"
binds=(--bind "$remote_v5/frozen/source:$repo:ro" --bind "$remote_v5:$local_v5:ro"
 --bind "$remote_addon:$local_addon:ro" --bind "$remote_addon/jobs:$local_addon/jobs:rw"
 --bind "$remote_addon/results:$local_addon/results:rw")
export NUMBA_CACHE_DIR="$visible_out/numba_cache" MPLCONFIGDIR="$visible_out/matplotlib"
run_bounded apptainer exec "${binds[@]}" --pwd "$local_addon" "$image" "$python" \
 "$local_v5/prepare_two_shock_torch.py" verify --stage "$local_v5" > "$out/container_verification.json"
run_bounded apptainer exec "${binds[@]}" --pwd "$local_addon" "$image" "$python" \
 "$local_addon/diagnose_initial_path.py" --config "$local_addon/config.json" --output "$visible_out/run" \
 > "$out/driver.log" 2>&1
