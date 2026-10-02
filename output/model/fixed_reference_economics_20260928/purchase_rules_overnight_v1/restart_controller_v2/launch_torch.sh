#!/usr/bin/env bash
#SBATCH --job-name=purchase80_restart
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/purchase_restart_controller_v2/logs/%x-%A_%a.out
set -euo pipefail
task=${SLURM_ARRAY_TASK_ID:?Specify reviewed stopped-chain IDs}
[[ "$task" =~ ^([0-9]|[1-3][0-9]|4[0-7])$ ]] || exit 2
restart=/scratch/td2248/projects/purchase_restart_controller_v2
cal=/scratch/td2248/projects/purchase_rules_overnight_v1
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
parent="$cal/results/chain_$task"
out="$restart/results/chain_$task"
[[ -f "$parent/launcher_start.json" && -f "$parent/search/search_completed.json" ]] || { echo 'Missing original receipts'; exit 2; }
"$python" - "$restart" "$cal" "$base" "$floor" <<'PY'
import hashlib,json,sys
from pathlib import Path
restart,cal,base,floor=map(Path,sys.argv[1:])
pins=json.loads((restart/'source_sha256.json').read_text())
for rel,digest in pins.items():
 p=restart/'source'/rel
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
for folder in (cal,base,floor):
 inv=json.loads((folder/'inventory.json').read_text())
 for rel,digest in inv['files'].items():
  p=folder/'source'/rel
  assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
PY
deadline_epoch=$("$python" - "$parent/launcher_start.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]));assert x['deadline_epoch']==x['start_epoch']+14400
print(x['deadline_epoch'])
PY
)
[[ $(date +%s) -lt $((deadline_epoch-900)) ]] || { echo 'No original search time remains'; exit 124; }
mkdir "$out"
mkdir "$out/numba_cache" "$out/matplotlib"
finish() {
 code=$?
 trap - EXIT
 "$python" - "$out" "$code" "$task" "$deadline_epoch" <<'PY'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),chain=int(sys.argv[3]),original_deadline_epoch=int(sys.argv[4]),finished_epoch=time.time()),indent=2)+'\n')
PY
 exit "$code"
}
trap finish EXIT
trap 'exit 143' TERM
export NUMBA_CACHE_DIR=/work/restart_results/numba_cache MPLCONFIGDIR=/work/restart_results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro" --bind "$parent:/work/parent:ro" --bind "$out:/work/restart_results:rw" --bind "$restart/source:/work/restart_source:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']): print(rel)
PY
)
binds+=(--bind "$cal/source/$packet:$repo/$packet:ro")
while IFS= read -r rel; do binds+=(--bind "$cal/source/$rel:$repo/$rel:ro"); done < <("$python" - "$cal/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if not rel.startswith(sys.argv[2]+'/'): print(rel)
PY
)
run_stage() {
 remaining=$((deadline_epoch-$(date +%s)-15))
 [[ "$remaining" -gt 0 ]] || exit 124
 timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" env RESTART_PACKET_ROOT="$repo/$packet" "$python" "$@"
}
run_stage /work/restart_source/controller.py --chain "$task" --parent /work/parent --out /work/restart_results > "$out/search.log" 2>&1
winner=$("$python" - "$out/restart_summary.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['winner'])
PY
)
if [[ "$winner" == restart ]]; then
 run_stage "$repo/$packet/run_psi.py" --chain "$task" --out /work/restart_results/postcheck --deadline-epoch "$deadline_epoch" --verify-only /work/restart_results/search/search_completed.json > "$out/postcheck.log" 2>&1
fi
