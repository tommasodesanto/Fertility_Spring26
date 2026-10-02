#!/usr/bin/env bash
#SBATCH --job-name=purchase_native_pilot3
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/purchase_native_integration_pilot_v3/logs/%x-%A_%a.out
set -euo pipefail
task=${SLURM_ARRAY_TASK_ID:?Array task required}
[[ "$task" == 0 || "$task" == 1 ]] || exit 2
arm=hard; [[ "$task" == 1 ]] && arm=quarter
pilot=/scratch/td2248/projects/purchase_native_integration_pilot_v3
mechanism=/scratch/td2248/projects/purchase_mechanism_v1
floor=/scratch/td2248/projects/normalized_floor_calibration_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
read -r chain root completed_sha < <("$python" - "$pilot" "$arm" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]);arm=sys.argv[2]
m=json.loads((p/'selection/manifest.json').read_text())
assert m['schema']=='provisional_native_integration_pilot_v1' and not m['final_calibration'] and not m['production_policy']
override=p/'source/selected_runtime.py'
assert hashlib.sha256(override.read_bytes()).hexdigest()==m['selected_runtime_override_sha256']
run_override=p/'source/run_case.py'
assert hashlib.sha256(run_override.read_bytes()).hexdigest()==m['run_case_override_sha256']
entry=m['arms'][arm]; selected=p/'selection'/f'selected_{arm}.json'
assert hashlib.sha256(selected.read_bytes()).hexdigest()==entry['selected_sha256']
row=json.loads(selected.read_text())
assert row['status']=='postchecked' and row['arm']==arm and row['chain']==entry['chain']
assert row['target_fingerprint']==m['target_fingerprint'] and row['weight_fingerprint']==m['weight_fingerprint']
assert row['remote_root']==entry['physical_root'] and row['origin']==entry['origin']
root=p/'selected_postchecks'/f'chain_{row["chain"]}' if entry['local_upload_required'] else Path(row['remote_root'])
assert root.is_dir()
completed=root/'postcheck/completed.json'
assert hashlib.sha256(completed.read_bytes()).hexdigest()==entry['completed_sha256']
assert json.loads(completed.read_text())['status']=='selected_numerically_verified'
report=root/'postcheck/selected_postcheck/phase_b_ge/selected_root'
for name,digest in row['report_sha256'].items():
 assert hashlib.sha256((report/name).read_bytes()).hexdigest()==digest,name
arrays=root/'postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz'
assert arrays.stat().st_size==row['native_arrays_bytes']
print(row['chain'],root,entry['completed_sha256'])
PY
)
"$python" - "$mechanism" "$base" "$floor" <<'PY'
import hashlib,json,sys
from pathlib import Path
for value in sys.argv[1:]:
 root=Path(value);pins=json.loads((root/'inventory.json').read_text())['files']
 for rel,digest in pins.items():
  assert hashlib.sha256((root/'source'/rel).read_bytes()).hexdigest()==digest,rel
PY
out="$pilot/results/$arm"
mkdir "$out" || { echo 'Refusing existing pilot result'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+1140))
terminal() {
 code=$?; trap - EXIT
 "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" "$chain" "$arm" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),start_epoch=int(sys.argv[3]),deadline_epoch=int(sys.argv[4]),chain=int(sys.argv[5]),arm=sys.argv[6],pilot=True,finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
 exit "$code"
}
trap terminal EXIT
trap 'exit 143' TERM
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']): print(rel)
PY
)
binds+=(--bind "$mechanism/source/$packet:$repo/$packet:ro")
binds+=(--bind "$pilot/source/selected_runtime.py:$repo/$packet/mechanism/selected_runtime.py:ro")
binds+=(--bind "$pilot/source/run_case.py:$repo/$packet/mechanism/run_case.py:ro")
binds+=(--bind "$mechanism/source/code/model/experiments/transition_readiness:$repo/code/model/experiments/transition_readiness:ro")
binds+=(--bind "$mechanism/source/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json:$repo/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json:ro")
binds+=(--bind "$root:$repo/$packet/results/chain_${chain}:ro")
binds+=(--bind "$pilot/selection:$repo/$packet/collection/readout:ro")
remaining=$((deadline_epoch-$(date +%s)-15))
(( remaining > 0 )) || exit 124
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/mechanism/run_case.py" \
 --arm "$arm" --selected-json "$repo/$packet/collection/readout/selected_${arm}.json" \
 --selected-completed "$repo/$packet/results/chain_${chain}/postcheck/completed.json" \
 --kind control --horizon 1 --smoke-one-date --out /work/results/run \
 --deadline-epoch "$deadline_epoch" --maximum-policy-calls 64 > "$out/run.log" 2>&1
