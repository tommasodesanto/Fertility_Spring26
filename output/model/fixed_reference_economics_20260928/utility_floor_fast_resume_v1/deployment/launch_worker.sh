#!/usr/bin/env bash
#SBATCH --job-name=utility_floor_fast_resume_v1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/utility_floor_fast_resume_v2/logs/%x-%A_%a.out
set -euo pipefail
start_epoch=$(date +%s)
remote=/scratch/td2248/projects/utility_floor_fast_resume_v2
old=/scratch/td2248/projects/utility_calibration_round1_v2
deadline_epoch=${FAST_DEADLINE_EPOCH:?}
lane=${FAST_STAGE:?}_${SLURM_ARRAY_TASK_ID:?}
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/utility_calibration_round1_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
out="$FAST_RESULT_ROOT/$FAST_STAGE/${SLURM_ARRAY_TASK_ID}_launcher"
mkdir "$out"
terminal_receipt() {
  local code=$?
  trap - EXIT
  "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" "$lane" <<'PY'
import json,os,sys,time
from pathlib import Path
out=Path(sys.argv[1]);code=int(sys.argv[2]);now=time.time()
r=dict(status='launcher_completed' if code==0 else 'launcher_failed',exit_code=code,lane=sys.argv[5],start_epoch=int(sys.argv[3]),deadline_epoch=float(sys.argv[4]),finished_epoch=now,elapsed_seconds=now-int(sys.argv[3]),slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'),no_automatic_restart=True)
p=out/'launcher_terminal.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(r,indent=2)+'\n');t.replace(p)
PY
  exit "$code"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
"$python" - "$old" "$base" "$inputs" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,base,inputs=map(Path,sys.argv[1:])
for folder in (base,remote):
    manifest=json.loads((folder/'inventory.json').read_text())
    for rel,digest in manifest['files'].items():
        assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest,rel
# Original runtime inventory is fully authenticated above.
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7','arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
    assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
PY
mkdir "$out/numba_cache" "$out/matplotlib"
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$remote/results:$remote/results:rw" --bind "$old/results/search/floor_s2/run:$old/results/search/floor_s2/run:ro" --bind "$remote/source:$repo/output/model/fixed_reference_economics_20260928/utility_floor_fast_resume_v1:ro" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$old/source/$packet:$repo/$packet:ro")
while IFS= read -r rel; do binds+=(--bind "$old/source/$rel:$repo/$rel:ro"); done < <("$python" - "$old/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
    if not rel.startswith(sys.argv[2]+'/'):print(rel)
PY
)
remaining=$("$python" - "$deadline_epoch" <<'PYD'
import sys,time
print(max(0,int(float(sys.argv[1])-time.time()-15)))
PYD
)
[[ "$remaining" -gt 0 ]] || exit 124
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/output/model/fixed_reference_economics_20260928/utility_floor_fast_resume_v1/worker.py" --runtime "$repo/$packet" --task-file "$FAST_TASK_FILE" --index "$SLURM_ARRAY_TASK_ID" --out "$FAST_RESULT_ROOT/$FAST_STAGE/$SLURM_ARRAY_TASK_ID" > "$out/native.log" 2>&1
