#!/usr/bin/env bash
#SBATCH --job-name=norm_hP_diag
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/normalized_floor_extension_v2/logs/%x-%A_%a.out
set -euo pipefail

mode=${FLOOR_MODE:?FLOOR_MODE must be gate or point}
[[ "$mode" == gate || "$mode" == point ]] || exit 2
index=${SLURM_ARRAY_TASK_ID:-0}
if [[ "$mode" == point ]]; then [[ "$index" =~ ^[0-3]$ ]] || exit 2; fi
remote=/scratch/td2248/projects/normalized_floor_extension_v2
v2=/scratch/td2248/projects/normalized_calibration_v2
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/normalized_floor_extension_v2
v2packet=output/model/fixed_reference_economics_20260928/normalized_calibration_v2
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
out="$remote/results/$mode$index"
mkdir "$out" || { echo 'Refusing existing results'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal() {
  code=$?
  trap - EXIT
  "$python" - "$out" "$code" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1], 'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$code"
}
trap terminal EXIT
trap 'exit 143' TERM
"$python" - "$out" "$remote/source/$packet/manifest.json" "$mode" "$index" <<'PY'
import hashlib,json,os,sys,time
from pathlib import Path
out,manifest,mode,index=sys.argv[1:]
Path(out,'launcher_start.json').write_text(json.dumps(dict(
    mode=mode,index=int(index),started_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),
    source_manifest_sha256=hashlib.sha256(Path(manifest).read_bytes()).hexdigest()),indent=2)+'\n')
PY
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$v2/source/$v2packet:$repo/$v2packet:ro")
while IFS= read -r rel; do binds+=(--bind "$v2/source/$rel:$repo/$rel:ro"); done < <("$python" - "$v2/inventory.json" "$v2packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if not rel.startswith(sys.argv[2]+'/'): print(rel)
PY
)
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
args=(--mode "$mode" --out /work/results/run)
if [[ "$mode" == point ]]; then
  gate="$remote/results/gate0/run/completed.json"
  "$python" - "$gate" "$remote/source/$packet/manifest.json" <<'PY'
import hashlib,json,sys
gate=json.load(open(sys.argv[1]));assert gate['status']=='incumbent_replay_passed'
assert gate['source_manifest_sha256']==hashlib.sha256(open(sys.argv[2],'rb').read()).hexdigest()
PY
  binds+=(--bind "$remote/results/gate0/run:/work/gate:ro")
  args+=(--index "$index")
  export FLOOR_GATE_RECEIPT=/work/gate/completed.json
fi
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
timeout --signal=TERM --kill-after=10s 1780s apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/run.py" "${args[@]}" > "$out/run.log" 2>&1
