#!/usr/bin/env bash
#SBATCH --job-name=purchase_regions
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:15:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/purchase_broader_regions_v1/logs/%x-%A_%a.out
set -euo pipefail
slot=${SLURM_ARRAY_TASK_ID:?Requires reviewed 0-15 array}
[[ "$slot" =~ ^([0-9]|1[0-5])$ ]] || exit 2
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+4500))
(( deadline_epoch<=1790932500 )) || deadline_epoch=1790932500
[[ "$deadline_epoch" -gt $((start_epoch+900)) ]] || { echo '05:15 New York cutoff leaves no search budget'; exit 124; }
region=/scratch/td2248/projects/purchase_broader_regions_v1
cal=/scratch/td2248/projects/purchase_rules_overnight_v1
restart=/scratch/td2248/projects/purchase_restart_controller_v2
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
"$python" - "$region" "$cal" "$base" "$floor" "$restart" <<'PY'
import hashlib,json,sys
from pathlib import Path
region,cal,base,floor,restart=map(Path,sys.argv[1:])
for rel,digest in json.loads((region/'source_sha256.json').read_text()).items():
 p=region/'source'/rel
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
assert hashlib.sha256((region/'launch_torch.sh').read_bytes()).hexdigest()==json.loads((region/'source_sha256.json').read_text())['launch_torch.sh']
review=json.loads((restart/'source_sha256.json').read_text())
assert hashlib.sha256((restart/'source/controller.py').read_bytes()).hexdigest()==review['controller.py']
for folder in (cal,base,floor):
 for rel,digest in json.loads((folder/'inventory.json').read_text())['files'].items():
  p=folder/'source'/rel
  assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
PY
if [[ "${REGIONS_PREFLIGHT:-0}" == 1 ]]; then out="$region/preflight/slot_$slot"; else out="$region/results/slot_$slot"; fi
mkdir "$out"
mkdir "$out/numba_cache" "$out/matplotlib"
finish() {
 code=$?
 trap - EXIT
 "$python" - "$out" "$code" "$slot" "$start_epoch" "$deadline_epoch" <<'PY'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),slot=int(sys.argv[3]),start_epoch=int(sys.argv[4]),deadline_epoch=int(sys.argv[5]),absolute_deadline_epoch=1790932500,finished_epoch=time.time(),maximum_objective_calls=80,final_reserve_seconds=900,no_auto_retry=True),indent=2)+'\n')
PY
 exit "$code"
}
trap finish EXIT
trap 'exit 143' TERM
export NUMBA_CACHE_DIR=/work/region_results/numba_cache MPLCONFIGDIR=/work/region_results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro" --bind "$out:/work/region_results:rw" --bind "$region/source:/work/regions_source:ro" --bind "$restart/source:/work/restart_source:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):print(rel)
PY
)
binds+=(--bind "$cal/source/$packet:$repo/$packet:ro")
while IFS= read -r rel; do binds+=(--bind "$cal/source/$rel:$repo/$rel:ro"); done < <("$python" - "$cal/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if not rel.startswith(sys.argv[2]+'/'):print(rel)
PY
)
run_stage() {
 remaining=$((deadline_epoch-$(date +%s)-15))
 [[ "$remaining" -gt 0 ]] || exit 124
 timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" env REGIONS_PACKET_ROOT="$repo/$packet" REGIONS_REVIEWED_CENSOR=/work/restart_source/controller.py "$python" /work/regions_source/driver.py --slot "$slot" --stage "$1" --root /work/region_results --deadline-epoch "$deadline_epoch"
}
if [[ "${REGIONS_PREFLIGHT:-0}" == 1 ]]; then
 run_stage init > "$out/init.log" 2>&1
 exit 0
fi
run_stage search > "$out/search.log" 2>&1
if "$python" - "$out/search/search_completed.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1])).get('selected') is not None else 1)
PY
then
 run_stage postcheck > "$out/postcheck.log" 2>&1
fi
