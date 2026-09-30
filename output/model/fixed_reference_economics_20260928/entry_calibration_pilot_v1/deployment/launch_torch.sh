#!/usr/bin/env bash
#SBATCH --job-name=entry_calibration_pilot_v1
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/entry_calibration_pilot_v1/logs/%x-%A_%a.out
set -euo pipefail
# Global budget starts at launcher entry and includes authentication/preflight/native smoke.
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+3600))
remote=/scratch/td2248/projects/entry_calibration_pilot_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
arms=(empirical_credit zero_wealth nonnegative_mean)
task=${SLURM_ARRAY_TASK_ID:?Requires array task 0, 1, or 2}
[[ "$task" =~ ^[012]$ ]] || { echo 'Invalid arm index'; exit 2; }
arm=${arms[$task]}
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
# mkdir is atomic: a rerun cannot reuse or overwrite an arm's cache/results.
mkdir -p "$remote/results"
out="$remote/results/$arm"
mkdir "$out" || { echo "Refusing existing arm results: $out"; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do
  binds+=(--bind "$base/source/$rel:$repo/$rel:ro")
done < "$base/mounts.txt"
# Bind the complete new packet, as in the validated comparison deployment.
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
# Optional compact dependencies outside this packet are mounted separately.
while IFS= read -r rel; do
  binds+=(--bind "$remote/source/$rel:$repo/$rel:ro")
done < <("$python" - "$remote/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
    if not rel.startswith(sys.argv[2]+'/'):
        print(rel)
PY
)
run_mode() {
  local mode=$1 output=$2 limit_epoch=$3 remaining
  remaining=$((limit_epoch-$(date +%s)))
  [[ "$remaining" -gt 0 ]] || { echo 'Global pilot deadline reached'; exit 124; }
  timeout --signal=KILL "${remaining}s" apptainer exec "${binds[@]}" \
    --pwd "$repo" "$image" "$python" "$repo/$packet/runner.py" \
    --mode "$mode" --arm "$arm" --out "/work/results/$output" \
    --deadline-seconds "$remaining" \
    > "$out/$output.log" 2>&1
}
preflight_end=$((start_epoch+300))
run_mode preflight preflight "$preflight_end"
# Runner must gate search on its native full-GE baseline and exact repeat.
# No retry/requeue path exists: any failed preflight or run ends this arm.
run_mode run run "$deadline_epoch"
