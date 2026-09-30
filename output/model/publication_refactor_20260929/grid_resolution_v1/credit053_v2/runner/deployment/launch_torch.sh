#!/usr/bin/env bash
#SBATCH --job-name=grid_resolution_credit053_v2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/grid_resolution_credit053_v2/logs/%x-%j.out
set -euo pipefail
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+2400))
remote=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
paired=/scratch/td2248/projects/small_credit_replication_v2/source/arms/indexed
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
[[ ! -e "$remote/results" ]] || { echo 'Refusing existing results'; exit 2; }
[[ -d "$frozen" && -d "$inputs" && -d "$paired" ]]
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
# Verify ALL compact copied sources before mounting. Verify original input externally.
"$python" - "$remote" "$inputs" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote=Path(sys.argv[1]); inputs=Path(sys.argv[2])
receipt=json.loads((remote/'inventory.json').read_text())
for rel,digest in receipt['files'].items():
    assert hashlib.sha256((remote/'source'/rel).read_bytes()).hexdigest()==digest,rel
assert hashlib.sha256((inputs/'bundle.json').read_bytes()).hexdigest()=='427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7'
PY
mkdir "$remote/results"
binds=(--bind "$frozen:$repo:ro" --bind "$remote/results:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
# Overlay read-only packages/pinned files; never write into frozen project.
while IFS= read -r rel; do
  binds+=(--bind "$remote/source/$rel:$repo/$rel:ro")
done < "$remote/mounts.txt"
run_mode() {
  mode=$1; output=$2; external_end=$3
  remaining=$((external_end-$(date +%s)))
  [[ "$remaining" -gt 0 ]]
  timeout --signal=KILL "${remaining}s" apptainer exec "${binds[@]}" \
    --pwd "$repo" "$image" "$python" "$repo/$packet/runner/run_comparison.py" "$mode" \
    --out "/work/results/$output" --deadline-epoch "$deadline_epoch" \
    > "$remote/results/$output.log" 2>&1
}
# Same subprocess/price-loop orchestration, mocked solve/observer; zero lifecycle calls.
run_mode preflight preflight "$((start_epoch+300))"
# Fresh distinct caches are created in each full arm output; no warm-cache reuse.
run_mode full full "$deadline_epoch"
"$python" - "$remote/results/full/completed.json" "$remote/results/preflight/completed.json" <<'PY'
import json,sys
full=json.load(open(sys.argv[1])); smoke=json.load(open(sys.argv[2]))
assert smoke['lifecycle_solves']==0 and smoke['status']=='mock_exact_loop_zero_solves'
assert full['status']=='full_passed' and full['lifecycle_solves']<=40
PY
end_epoch=$(date +%s)
"$python" - "$remote/results/job_workflow_timing.json" "$start_epoch" "$end_epoch" <<'PY'
import json,sys
json.dump(dict(start_epoch=int(sys.argv[2]),end_epoch=int(sys.argv[3]),complete_job_workflow_seconds=int(sys.argv[3])-int(sys.argv[2]),budget_seconds=2400,includes='verification, preflight, production imports/compilation/authentication/solves/observers/plots/comparison'),open(sys.argv[1],'w'),indent=2)
PY
