#!/usr/bin/env bash
#SBATCH --job-name=entry_ratio_comparison_v1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/entry_ratio_comparison_v1/logs/%x-%j.out
set -euo pipefail
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+2400))
remote=/scratch/td2248/projects/entry_ratio_comparison_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/entry_ratio_comparison_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
[[ ! -e "$remote/results" ]] || { echo 'Refusing existing results'; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
"$python" - "$remote" "$base" "$inputs" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,base,inputs=map(Path,sys.argv[1:])
for folder in (base,remote):
    manifest=json.loads((folder/'inventory.json').read_text())
    for rel,digest in manifest['files'].items():
        assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest,rel
assert hashlib.sha256((inputs/'bundle.json').read_bytes()).hexdigest()=='427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7'
PY
mkdir "$remote/results"
binds=(--bind "$frozen:$repo:ro" --bind "$remote/results:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do
  binds+=(--bind "$base/source/$rel:$repo/$rel:ro")
done < "$base/mounts.txt"
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
run_mode() {
  mode=$1; output=$2; external_end=$3
  remaining=$((external_end-$(date +%s)))
  [[ "$remaining" -gt 0 ]]
  timeout --signal=KILL "${remaining}s" apptainer exec "${binds[@]}" \
    --pwd "$repo" "$image" "$python" "$repo/$packet/run_comparison.py" "$mode" \
    --out "/work/results/$output" --deadline-epoch "$deadline_epoch" \
    > "$remote/results/$output.log" 2>&1
}
run_mode preflight preflight "$((start_epoch+300))"
run_mode full full "$deadline_epoch"
"$python" - "$remote/results/full/completed.json" "$remote/results/preflight/completed.json" <<'PY'
import json,sys
full=json.load(open(sys.argv[1]));smoke=json.load(open(sys.argv[2]))
assert full['status']=='full_passed' and full['lifecycle_solves']<=40
assert smoke['status']=='mock_exact_loop_zero_solves' and smoke['lifecycle_solves']==0
PY
