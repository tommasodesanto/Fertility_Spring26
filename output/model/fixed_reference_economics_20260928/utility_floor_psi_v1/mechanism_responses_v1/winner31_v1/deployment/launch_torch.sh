#!/usr/bin/env bash
#SBATCH --job-name=floor_mech_v1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/utility_floor_winner31_responses_v1/logs/%x-%j.out
set -euo pipefail
remote=/scratch/td2248/projects/utility_floor_winner31_responses_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+1200))
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
mode=${1:-run}
[[ "$mode" == run || "$mode" == --initializer-only ]] || exit 2
out="$remote/results"
[[ "$mode" != --initializer-only ]] || out="$remote/zeroLC_preflight"
mkdir "$out" || { echo 'Refusing existing results'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal() {
 local code=$?
 trap - EXIT
 "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" <<'PY'
import json,os,sys,time
from pathlib import Path
r=dict(exit_code=int(sys.argv[2]),start_epoch=int(sys.argv[3]),deadline_epoch=int(sys.argv[4]),finished_epoch=time.time(),launcher_pid=os.getppid(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True)
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(r,indent=2)+'\n')
PY
 exit "$code"
}
trap terminal EXIT
trap 'exit 143' TERM
"$python" - "$remote" "$base" "$inputs" "$out" "$start_epoch" "$deadline_epoch" "$$" <<'PY'
import hashlib,json,os,sys,time
from pathlib import Path
remote,base,inputs,out=map(Path,sys.argv[1:5]);manifest=json.loads((remote/'inventory.json').read_text())
for folder in (base,remote):
 for rel,digest in json.loads((folder/'inventory.json').read_text())['files'].items():
  assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest,rel
for rel,item in manifest['reused_large_inputs'].items():assert hashlib.sha256(Path(item['remote_path']).read_bytes()).hexdigest()==item['sha256'],rel
rel='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/mechanism_responses_v1/winner31_v1/deployment/launch_torch.sh'
assert hashlib.sha256((remote/'launch_torch.sh').read_bytes()).hexdigest()==manifest['files'][rel]
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7','arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
r=dict(status='source_hashes_passed',overlay_sources=len(manifest['files']),reused_large_inputs=manifest['reused_large_inputs'],start_epoch=int(sys.argv[5]),deadline_epoch=int(sys.argv[6]),launcher_pid=int(sys.argv[7]),slurm_job_id=os.getenv('SLURM_JOB_ID'),cpu=1,memory_gib=24,threads=1,wall_seconds=1200,no_auto_retry=True)
(out/'launcher_start.json').write_text(json.dumps(r,indent=2)+'\n')
PY
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
# Bind packet directories as a whole so every frozen local candidate path exists.
binds+=(--bind "$remote/source/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1:$repo/output/model/fixed_reference_economics_20260928/utility_floor_psi_v1:ro")
while IFS= read -r rel; do binds+=(--bind "$remote/source/$rel:$repo/$rel:ro"); done < "$remote/mounts.txt"
run_stage() {
 local name=$1 script=$2 remaining
 shift 2
 remaining=$((deadline_epoch-$(date +%s)-15))
 [[ "$remaining" -gt 0 ]] || { echo 'Absolute20min deadline reached'; exit 124; }
 timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$script" "$@" > "$out/$name.log" 2>&1
}
# Source/input and actual runtime gates are inside the same absolute budget; zero extra native smoke solves.
run_stage initializer "$packet/deployment/initialize_runtime.py" --out /work/results/initializer
[[ "$mode" != --initializer-only ]] || exit 0
run_stage responses "$packet/fixed_price_responses.py" --out /work/results/responses --deadline-epoch "$deadline_epoch"
