#!/usr/bin/env bash
#SBATCH --job-name=estate_parallel_candidate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/current_estate_parallel_candidates_20261004_v1/logs/%x-%j.out
set -euo pipefail

# One invocation owns one candidate. The add-on must have been staged and pinned
# before submission; neither this launcher nor Slurm creates a replacement run.
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
local_addon="$repo/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/parallel_candidates_v1"
local_v5="$repo/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5"
remote_v5=/scratch/td2248/projects/current_estate_two_shock_20261004_v5
expected_remote_addon=/scratch/td2248/projects/current_estate_parallel_candidates_20261004_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif

[[ "$#" == 3 ]] || { echo 'Usage: launch_parallel_transition_candidate.sh REMOTE_ADDON TASK_ID CONFIG_SHA256' >&2; exit 2; }
remote_addon=$1
task_id=$2
expected_config_sha=$3
[[ "$remote_addon" == "$expected_remote_addon" ]] || { echo 'Unrecognized remote add-on root' >&2; exit 2; }
[[ "$task_id" =~ ^c0[0-5]_h(24|32)$ ]] || { echo 'Unrecognized candidate task ID' >&2; exit 2; }
[[ "$expected_config_sha" =~ ^[[:xdigit:]]{64}$ ]] || { echo 'Expected exact 64-character config SHA256' >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" && -z "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'A non-array Slurm job is required' >&2; exit 2; }
[[ "${SLURM_NTASKS:-1}" == 1 && "${SLURM_CPUS_PER_TASK:-1}" == 1 ]] || { echo 'Exactly one task and one CPU required' >&2; exit 2; }
[[ "${SLURM_MEM_PER_NODE:-24576}" =~ ^[0-9]+$ && "${SLURM_MEM_PER_NODE:-24576}" -le 24576 ]] || { echo '24 GiB memory cap exceeded or unknown' >&2; exit 2; }
[[ -d "$remote_addon" && -d "$remote_addon/logs" && -d "$remote_addon/jobs" && -d "$remote_addon/results" ]] || { echo 'Pre-staged add-on root is incomplete' >&2; exit 2; }
[[ -d "$remote_v5" ]] || { echo 'Immutable v5 runtime is absent' >&2; exit 2; }

module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH

config="$remote_addon/configs/$task_id.json"
driver="$remote_addon/parallel_transition_candidate.py"
[[ -f "$config" && -f "$driver" ]] || { echo 'Candidate config or driver is absent' >&2; exit 2; }

# mkdir is the atomic no-retry claim. A failed candidate retains its claim.
claim="$remote_addon/jobs/$task_id.claim"
mkdir "$claim" || { echo "Refusing existing or unknown candidate claim: $claim" >&2; exit 2; }
printf '%s\n' "$SLURM_JOB_ID" > "$claim/slurm_job_id"
printf '%s\n' "$expected_config_sha" > "$claim/config_sha256"
job="$remote_addon/jobs/${task_id}_${SLURM_JOB_ID}"
out="$remote_addon/results/$task_id"
mkdir "$job" || { echo "Refusing existing or unknown job: $job" >&2; exit 2; }
mkdir "$out" || { echo "Refusing existing or unknown candidate output: $out" >&2; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
visible_out="$local_addon/results/$task_id"
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+10800))
child_pid=''
heartbeat_pid=''

stop_child() {
 if [[ -n "$child_pid" ]] && kill -0 "$child_pid" 2>/dev/null; then
  kill -TERM -- "-$child_pid" 2>/dev/null || true
  for _ in {1..10}; do kill -0 "$child_pid" 2>/dev/null || break; sleep 1; done
  kill -KILL -- "-$child_pid" 2>/dev/null || true
  wait "$child_pid" 2>/dev/null || true
 fi
}
terminal() {
 local status=$?
 trap - EXIT TERM INT
 stop_child
 if [[ -n "$heartbeat_pid" ]]; then
  kill -TERM -- "-$heartbeat_pid" 2>/dev/null || true
  wait "$heartbeat_pid" 2>/dev/null || true
 fi
 "$python" - "$job/launcher_terminal.json" "$status" "$start_epoch" "$deadline_epoch" <<'PYRECEIPT'
import json, os, sys, time
from pathlib import Path
path, status, start, deadline = sys.argv[1:]
Path(path).write_text(json.dumps(dict(status='complete' if int(status)==0 else 'failed',
    exit_code=int(status), start_epoch=int(start), deadline_epoch=int(deadline),
    finished_epoch=time.time(), slurm_job_id=os.getenv('SLURM_JOB_ID'),
    node=os.getenv('SLURMD_NODENAME'), no_auto_retry=True), sort_keys=True, indent=2)+'\n')
PYRECEIPT
 cp "$job/launcher_terminal.json" "$out/launcher_terminal.json"
 exit "$status"
}
trap terminal EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

# The submitted SHA authenticates the whole config. Its driver pin authenticates
# the separately staged add-on program before even constructing the model.
"$python" - "$config" "$driver" "$local_addon/parallel_transition_candidate.py" "$expected_config_sha" "$task_id" > "$out/addon_verification.json" <<'PYVERIFY'
import hashlib, json, sys
from pathlib import Path
config_path, driver_path, expected_driver_path, expected_config, task_id = sys.argv[1:]
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
actual_config = sha(config_path)
if actual_config != expected_config:
    raise SystemExit('Remote config SHA256 differs from submission pin')
config = json.loads(Path(config_path).read_text())
pin = config.get('driver')
if not isinstance(pin, dict) or set(pin) != {'path', 'sha256'} or pin['path'] != expected_driver_path or pin['sha256'] != sha(driver_path):
    raise SystemExit('Driver bytes differ from the exact config driver pin')
if config.get('task_id') != task_id:
    raise SystemExit('Config task_id differs from submitted task ID')
print(json.dumps(dict(task_id=task_id, config_sha256=actual_config,
    driver_sha256=sha(driver_path), verified=True), sort_keys=True))
PYVERIFY

"$python" - "$job/launcher_start.json" "$start_epoch" "$deadline_epoch" "$expected_config_sha" "$task_id" <<'PYRECEIPT'
import json, os, sys
from pathlib import Path
path, start, deadline, sha, task_id = sys.argv[1:]
Path(path).write_text(json.dumps(dict(task_id=task_id, start_epoch=int(start),
    deadline_epoch=int(deadline), config_sha256=sha,
    output=f'results/{task_id}/run', cpus=1, memory_gib=24,
    numba_threads=1, blas_threads=1, openmp_threads=1,
    slurm_job_id=os.getenv('SLURM_JOB_ID'), node=os.getenv('SLURMD_NODENAME'),
    no_auto_retry=True), sort_keys=True, indent=2)+'\n')
PYRECEIPT
cp "$job/launcher_start.json" "$out/launcher_start.json"
setsid bash -c 'while true; do date -u +%Y-%m-%dT%H:%M:%SZ > "$1"; sleep 60; done' \
 _ "$job/launcher_heartbeat.txt" >/dev/null 2>&1 < /dev/null &
heartbeat_pid=$!
run_bounded() {
 local remaining=$((deadline_epoch-$(date +%s)-60))
 [[ "$remaining" -gt 0 ]] || return 124
 setsid timeout --signal=TERM --kill-after=10s "${remaining}s" "$@" <&0 & child_pid=$!
 wait "$child_pid"
 local status=$?
 child_pid=''
 return "$status"
}

# The frozen v5 verifier checks the full host inventory (4,599 files); repeat
# inside the container to authenticate the actual read-only source mount.
run_bounded "$python" "$remote_v5/prepare_two_shock_torch.py" verify --stage "$remote_v5" > "$out/host_verification.json"
binds=(--bind "$remote_v5/frozen/source:$repo:ro" --bind "$remote_v5:$local_v5:ro"
 --bind "$remote_addon:$local_addon:ro" --bind "$remote_addon/jobs:$local_addon/jobs:rw"
 --bind "$remote_addon/results:$local_addon/results:rw")
export NUMBA_CACHE_DIR="$visible_out/numba_cache" MPLCONFIGDIR="$visible_out/matplotlib"
run_bounded apptainer exec "${binds[@]}" --pwd "$local_addon" "$image" "$python" \
 "$local_v5/prepare_two_shock_torch.py" verify --stage "$local_v5" > "$out/container_verification.json"
run_bounded apptainer exec "${binds[@]}" --pwd "$local_addon" "$image" "$python" - \
 "$local_addon/configs/$task_id.json" "$local_addon/parallel_transition_candidate.py" "$expected_config_sha" <<'PYCONTAINER' > "$out/container_addon_verification.json"
import hashlib,json,sys
from pathlib import Path
c,d,pin=sys.argv[1:]
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
cfg=json.loads(Path(c).read_text())
if sha(c)!=pin or cfg.get('driver')!={'path':d,'sha256':sha(d)}:
    raise SystemExit('Container add-on authentication failed')
print(json.dumps(dict(config_sha256=sha(c),driver_sha256=sha(d),verified=True),sort_keys=True))
PYCONTAINER
run_bounded apptainer exec "${binds[@]}" --pwd "$local_addon" "$image" "$python" \
 "$local_addon/parallel_transition_candidate.py" --config "$local_addon/configs/$task_id.json" \
 --output "$visible_out/run" > "$out/driver.log" 2>&1
