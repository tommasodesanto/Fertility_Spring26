#!/usr/bin/env bash
# Immutable source overlay. Invoke "preflight" outside Slurm for zero-solve checks.
# Slurm modes: ESTATE_RUN_MODE=smoke|production. Production is submitted separately.
#SBATCH --job-name=estatebirth
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/estate_birth_overnight_20261004_v1/logs/%x-%A_%a.out
set -euo pipefail

remote=/scratch/td2248/projects/estate_birth_overnight_20261004_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
floor_remote=/scratch/td2248/projects/normalized_floor_calibration_v1
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
mode=${1:-${ESTATE_RUN_MODE:-production}}
[[ "$mode" =~ ^(preflight|mock|production)$ ]] || { echo "Invalid mode: $mode"; exit 2; }
if [[ "$mode" == preflight ]]; then
  task=${PREFLIGHT_TASK:-0}
else
  task=${SLURM_ARRAY_TASK_ID:?Requires Slurm array task 0..19}
fi
[[ "$task" =~ ^(0|[1-9]|1[01])$ ]] || { echo 'Invalid continuation chain'; exit 2; }
arm=binary; chain=$task
starts_file=/work/deployment/control/starts.json
starts_sha=$("$python" - "$remote/control/starts.json" <<'PY'
import hashlib,sys
print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())
PY
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
if [[ "$mode" == preflight ]]; then
  apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/verify_stage.py --container
  preflight_out="$remote/preflight/${arm}_chain_${chain}"
  mkdir -p "$preflight_out"
  timeout --signal=TERM --kill-after=10s 600s apptainer exec "${binds[@]}" \
    --bind "$preflight_out:/work/results:rw" --pwd "$repo" "$image" "$python" \
    "$repo/code/model/experiments/birth_count_choice/cluster_calibrate.py" \
    --arm "$arm" --chain "$chain" --out /work/results/run \
    --deadline-epoch "$(($(date +%s)+3600))" --starts-file "$starts_file" \
    --starts-file-sha256 "$starts_sha" --preflight-evaluator \
    > "$preflight_out/driver.log" 2>&1
  "$python" - "$preflight_out/run/completed.json" <<'PY'
import json,sys
assert json.load(open(sys.argv[1]))['status']=='evaluator_initialized_zero_solves'
PY
  exit 0
fi

start_epoch=$(date +%s)
if [[ "$mode" == production ]]; then wall_seconds=43200; else wall_seconds=300; fi
deadline_epoch=$((start_epoch+wall_seconds))
mkdir -p "$remote/results" "$remote/logs"
out="$remote/results/${mode}_binary_chain_${chain}"
mkdir "$out" || { echo "Refusing existing result: $out"; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal_receipt() {
  local status=$?
  trap - EXIT
  "$python" - "$out" "$status" "$mode" "$arm" "$chain" "$start_epoch" "$deadline_epoch" <<'PY'
import json,os,sys,time
from pathlib import Path
out,status,mode,arm,chain,start,deadline=sys.argv[1:]
Path(out,'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(status),mode=mode,arm=arm,
    chain=int(chain),start_epoch=int(start),deadline_epoch=int(deadline),finished_epoch=time.time(),
    slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$status"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM
"$python" - "$out" "$mode" "$arm" "$chain" "$start_epoch" "$deadline_epoch" "$wall_seconds" "$remote" "$starts_sha" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
out,mode,arm,chain,start,deadline,wall,remote,starts_sha=sys.argv[1:]
stage_inventory_sha256=hashlib.sha256((Path(remote)/'inventory.json').read_bytes()).hexdigest()
Path(out,'launcher_start.json').write_text(json.dumps(dict(mode=mode,arm=arm,chain=int(chain),
    start_epoch=int(start),deadline_epoch=int(deadline),slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),
    cpus=1,memory_GiB=24,wall_seconds=int(wall),maximum_objective_calls=2 if mode=='smoke' else 500,starts_sha256=starts_sha,
    final_native_reserve_seconds=1800,stage_inventory_sha256=stage_inventory_sha256),indent=2)+'\n')
PY
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds+=(--bind "$out:/work/results:rw")
apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/verify_stage.py --container
options=(--starts-file "$starts_file" --starts-file-sha256 "$starts_sha")
[[ "$mode" == mock ]] && options+=(--mock-smoke)
remaining=$((deadline_epoch-$(date +%s)-15))
[[ "$remaining" -gt 0 ]] || exit 124
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
  "$repo/code/model/experiments/birth_count_choice/cluster_calibrate.py" \
  --arm "$arm" --chain "$chain" --out /work/results/run --deadline-epoch "$deadline_epoch" "${options[@]}" \
  > "$out/driver.log" 2>&1

if [[ "$mode" == smoke ]]; then
  "$python" - "$out/run/completed.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d['status']=='selected_numerically_verified' and d['objective_calls']==2
assert d['selected_postcheck']['status']=='passed' and d['repeat']['status']=='exact_full_ge_repeat_passed'
PY
fi
