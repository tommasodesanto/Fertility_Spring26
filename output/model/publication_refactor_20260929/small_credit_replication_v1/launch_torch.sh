#!/usr/bin/env bash
#SBATCH --job-name=small_credit_pair_v1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/small_credit_replication_v2/logs/%x-%j.out
set -euo pipefail
start_epoch="$(date +%s)"
deadline_epoch="$((start_epoch+2400))"
remote=/scratch/td2248/projects/small_credit_replication_v2
source="$remote/source"
results="$remote/results"
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
[[ ! -e "$results/scalar" && ! -e "$results/indexed" ]] || { echo 'Refusing duplicate'; exit 2; }
(cd "$source" && sha256sum -c source.sha256)
[[ "$(sha256sum "$inputs/bundle.json" | cut -d' ' -f1)" == 427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7 ]]
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
mkdir -p "$results" "$remote/cache/scalar" "$remote/cache/indexed"
run_arm() {
  arm="$1"; mode="$2"; out="$3"; until="$4"
  export NUMBA_CACHE_DIR=/work/cache
  export PYTHONPATH=/work/source/arms/$arm:/work/source/arms/$arm/source
  remaining="$((until-$(date +%s)))"
  [[ "$remaining" -gt 0 ]]
  timeout --signal=TERM --kill-after=5s "${remaining}s" apptainer exec \
    --bind "$frozen:$repo:ro,$source:/work/source:ro,$inputs:/work/inputs:ro,$results:/work/results:rw,$remote/cache/$arm:/work/cache:rw" \
    --pwd "$repo" /share/apps/images/ubuntu-24.04.4.sif \
    /share/apps/anaconda3/2025.06/bin/python /work/source/arms/$arm/run_arm.py "$mode" \
    --reference-root "$repo" --bundle /work/inputs --out "/work/results/$out" --deadline-epoch "$deadline_epoch" \
    > "$results/$out.log" 2>&1
}
# Authenticate and traverse both loops with zero lifecycle evaluations first.
for arm in scalar indexed; do run_arm "$arm" smoke "preflight_$arm" "$((start_epoch+300))"; done
# Preflight imports can populate Numba caches. Production arms receive fresh directories.
mkdir "$remote/cache/scalar_production" "$remote/cache/indexed_production"
# Rename the preflight caches; replace arm paths with their empty production directories.
mv "$remote/cache/scalar" "$remote/cache/scalar_preflight"
mv "$remote/cache/indexed" "$remote/cache/indexed_preflight"
mv "$remote/cache/scalar_production" "$remote/cache/scalar"
mv "$remote/cache/indexed_production" "$remote/cache/indexed"
for arm in scalar indexed; do
  arm_start="$(date +%s)"
  run_arm "$arm" full "$arm" "$deadline_epoch"
  arm_end="$(date +%s)"
  printf '{"elapsed_seconds":%s}\n' "$((arm_end-arm_start))" > "$results/$arm/workflow_timing.json"
  /share/apps/anaconda3/2025.06/bin/python - "$results/$arm/completed.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1])); assert r['status']=='full_passed' and r['lifecycle_solves']==6
PY
done
export PYTHONPATH=
/share/apps/anaconda3/2025.06/bin/python "$source/compare_pair.py" "$results/scalar" "$results/indexed" "$results/comparison.json" /scratch/td2248/projects/small_credit_v1/results/full
