#!/usr/bin/env bash
#SBATCH --job-name=parenthood_floor_quick_v1
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/parenthood_floor_quick_v1/logs/%x-%A_%a.out
set -euo pipefail
start_epoch=$(date +%s)
common_deadline_epoch=1790811660 # 2026-09-30T23:41:00Z
deadline_epoch=$((start_epoch+1800))
if [[ "$common_deadline_epoch" -lt "$deadline_epoch" ]]; then deadline_epoch=$common_deadline_epoch; fi
remote=/scratch/td2248/projects/parenthood_floor_quick_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
lanes=(floor no_A constant_alpha)
task=${SLURM_ARRAY_TASK_ID:?Requires array task 0, 1, or 2}
[[ "$task" =~ ^[012]$ ]] || { echo 'Invalid arm index'; exit 2; }
lane=${lanes[$task]}
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
mkdir -p "$remote/results"
if [[ "${CAL_PREFLIGHT_ONLY:-0}" == 1 ]]; then
  mkdir -p "$remote/preflight_validation"
  out="$remote/preflight_validation/$lane"
else
  out="$remote/results/$lane"
fi
mkdir "$out" || { echo "Refusing existing lane results: $out"; exit 2; }
terminal_receipt() {
  local code=$?
  trap - EXIT
  "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" "$lane" <<'PY'
import json,os,sys,time
from pathlib import Path
out=Path(sys.argv[1]);code=int(sys.argv[2]);now=time.time()
r=dict(status='launcher_completed' if code==0 else 'launcher_failed',exit_code=code,lane=sys.argv[5],start_epoch=int(sys.argv[3]),deadline_epoch=int(sys.argv[4]),finished_epoch=now,elapsed_seconds=now-int(sys.argv[3]),slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_array_job_id=os.getenv('SLURM_ARRAY_JOB_ID'),slurm_array_task_id=os.getenv('SLURM_ARRAY_TASK_ID'),no_automatic_restart=True)
p=out/'launcher_terminal.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(r,indent=2)+'\n');t.replace(p)
PY
  exit "$code"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
"$python" - "$remote" "$base" "$inputs" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,base,inputs=map(Path,sys.argv[1:])
for folder in (base,remote):
    manifest=json.loads((folder/'inventory.json').read_text())
    for rel,digest in manifest['files'].items():
        assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest,rel
launcher_rel='output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/deployment/launch_torch.sh'
assert hashlib.sha256((remote/'launch_torch.sh').read_bytes()).hexdigest()==json.loads((remote/'inventory.json').read_text())['files'][launcher_rel], 'Entrypoint launcher drift'
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7','arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
    assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
PY
mkdir "$out/numba_cache" "$out/matplotlib"
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
while IFS= read -r rel; do binds+=(--bind "$remote/source/$rel:$repo/$rel:ro"); done < <("$python" - "$remote/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
    if not rel.startswith(sys.argv[2]+'/'):print(rel)
PY
)
run_mode() {
  local mode=$1 output=$2 limit_epoch=$3 remaining
  remaining=$((limit_epoch-$(date +%s)-15))
  [[ "$remaining" -gt 0 ]] || { echo 'Global calibration deadline reached'; exit 124; }
  timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" \
    --pwd "$repo" "$image" "$python" "$repo/$packet/runner.py" \
    --mode "$mode" --arm "$lane" --out "/work/results/$output" \
    --deadline-seconds "$remaining" \
    > "$out/$output.log" 2>&1
}
if [[ "${CAL_PREFLIGHT_ONLY:-0}" == 1 ]]; then
  preflight_end=$(($(date +%s)+300))
  if [[ "$deadline_epoch" -lt "$preflight_end" ]]; then preflight_end=$deadline_epoch; fi
  run_mode preflight preflight "$preflight_end"
  exit 0
fi
preflight_end=$((start_epoch+300))
if [[ "$deadline_epoch" -lt "$preflight_end" ]]; then preflight_end=$deadline_epoch; fi
run_mode preflight preflight "$preflight_end"
# Native experiment and reporting share the common 30-minute clock.
run_mode run run "$deadline_epoch"
