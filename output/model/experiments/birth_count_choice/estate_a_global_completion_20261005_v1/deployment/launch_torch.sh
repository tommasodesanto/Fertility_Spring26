#!/usr/bin/env bash
#SBATCH --job-name=estateglobal
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/estate_birth_global_completion_20261005_v1/logs/%x-%A_%a.out
set -euo pipefail
remote=/scratch/td2248/projects/estate_birth_global_completion_20261005_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
floor_remote=/scratch/td2248/projects/normalized_floor_calibration_v1
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
mode=${ESTATE_RUN_MODE:?Set ESTATE_RUN_MODE=smoke|preflight|production}
[[ "$mode" =~ ^(smoke|preflight|production)$ ]] || exit 2
task=${SLURM_ARRAY_TASK_ID:?Requires Slurm task}
[[ "$task" =~ ^(0|[1-8])$ ]] || exit 2
plan_sha=$("$python" - "$remote/control/plan.json" <<'PY2'
import hashlib,sys
print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())
PY2
)
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
"$python" "$remote/verify_stage.py" --host
"$python" - "$inputs" <<'PY'
import hashlib,sys
from pathlib import Path
root=Path(sys.argv[1])
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7',
                    'arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
PY
for source_root in "$base" "$floor_remote"; do
  "$python" - "$source_root" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]); inv=json.loads((root/'inventory.json').read_text())
for rel,digest in inv['files'].items():
    path=root/'source'/rel
    assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==digest,rel
PY
done

parent_remote=/scratch/td2248/projects/estate_birth_calibration_20261003_v3
binds=(--bind "$frozen:$repo:ro" --bind "$remote:/work/deployment:ro" --bind "/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2:/work/parent_v2:ro" --bind "$parent_remote:/work/parent:ro"
       --bind "/scratch/td2248/projects/estate_birth_recovery_20261004_v1:/work/recovery:ro"
       --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor_remote/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor_remote/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']): print(rel)
PY
)
# Directory mounts give the new sandbox and result packets stable absolute paths.
for rel in \
  code/model/experiments/purchase_timing_sandbox \
  code/model/experiments/birth_count_choice \
  code/model/production/reference_inputs \
  output/model/fixed_reference_economics_20260928/purchase_timing_sandbox_v1 \
  output/model/fixed_reference_economics_20260928/soft_timing_review_v1 \
  output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1 \
  output/model/experiments/birth_count_choice/estate_a_calibration_v1 \
  output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1; do
  binds+=(--bind "$remote/source/$rel:$repo/$rel:ro")
done
# Apply all source pins last, covering any original files changed after the frozen base.
while IFS= read -r rel; do binds+=(--bind "$remote/source/$rel:$repo/$rel:ro"); done < <("$python" - "$remote/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if rel.startswith(('code/model/experiments/purchase_timing_sandbox/',
                    'code/model/experiments/birth_count_choice/',
                    'code/model/production/reference_inputs/',
                    'output/model/fixed_reference_economics_20260928/purchase_timing_sandbox_v1/',
                    'output/model/fixed_reference_economics_20260928/soft_timing_review_v1/',
                    'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/',
                    'output/model/experiments/birth_count_choice/estate_a_calibration_v1/',
                    'output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/')): continue
 print(rel)
PY
)
start_epoch=$(date +%s)
if [[ "$mode" == production ]]; then wall_seconds=5400; elif [[ "$mode" == preflight ]]; then wall_seconds=1500; else wall_seconds=300; fi
deadline_epoch=$((start_epoch+wall_seconds))
if [[ "$mode" == production && "$deadline_epoch" -gt 1791223200 ]]; then deadline_epoch=1791223200; fi
out="$remote/results/${mode}_task_${task}"
mkdir -p "$remote/results" "$remote/logs"
mkdir "$out" || { echo "Refusing existing result: $out"; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal_receipt() {
  local status=$?
  trap - EXIT
  "$python" - "$out" "$status" "$mode" "$task" <<'PY2'
import json,os,sys,time
from pathlib import Path
out,status,mode,task=sys.argv[1:]
Path(out,'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(status),mode=mode,task=int(task),finished_epoch=time.time(),
  slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY2
  exit "$status"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM
"$python" - "$out" "$mode" "$task" "$start_epoch" "$deadline_epoch" "$plan_sha" "$remote" <<'PY2'
import hashlib,json,os,sys
from pathlib import Path
out,mode,task,start,deadline,plan_sha,remote=sys.argv[1:]
Path(out,'launcher_start.json').write_text(json.dumps(dict(mode=mode,task=int(task),start_epoch=int(start),deadline_epoch=int(deadline),
  plan_sha256=plan_sha,stage_manifest_sha256=hashlib.sha256((Path(remote)/'stage_manifest.json').read_bytes()).hexdigest(),
  slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),cpus=1,memory_GiB=24),indent=2)+'\n')
PY2
# Keep recovery successors from overlapping the dependent global stage.
if [[ "$mode" == production ]]; then
  while true; do
    queue_snapshot=$(squeue -h -u "$USER" -o '%j %T') || { echo 'Recovery queue query failed; refusing native cases' >&2; exit 4; }
    if ! awk '$1 ~ /^(estatebirth|softtiming)$/ && $2 ~ /^(PENDING|RUNNING|COMPLETING)$/ {found=1} END {exit !found}' <<< "$queue_snapshot"; then break; fi
    "$python" - "$out" <<'PYWAIT'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],'queue_wait.json').write_text(json.dumps(dict(status='waiting_for_overnight_recovery_successors',epoch=time.time()))+'\n')
PYWAIT
    if (( $(date +%s) >= deadline_epoch - 300 )); then
      "$python" - "$out" <<'PYDEFER'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],'deferred_concurrency.json').write_text(json.dumps(dict(status='deferred_no_native_cases',epoch=time.time()))+'\n')
PYDEFER
      exit 0
    fi
    sleep 60
  done
fi
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds+=(--bind "$out:/work/results:rw")
apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/verify_stage.py --container
remaining=$((deadline_epoch-$(date +%s)-15))
[[ "$remaining" -gt 0 ]] || exit 124
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python"   /work/deployment/explore.py --mode "$mode" --task "$task" --out /work/results/run   --deadline-epoch "$deadline_epoch" --plan-sha256 "$plan_sha" > "$out/driver.log" 2>&1
"$python" - "$out/run/completed.json" "$mode" <<'PY2'
import json,sys
x=json.load(open(sys.argv[1]));mode=sys.argv[2]
if mode=='smoke': assert x['status']=='smoke_loop_passed_zero_solves' and x['completed_cases']==2
if mode=='preflight': assert x['status']=='native_preflight_passed' and x['completed_cases']==1
PY2
