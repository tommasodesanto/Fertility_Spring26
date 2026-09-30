#!/usr/bin/env bash
#SBATCH --job-name=single_market_verify_v1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodelist=cs713
#SBATCH --output=/scratch/td2248/projects/single_market_verification_v1/logs/%x-%j.out
set -euo pipefail
start_epoch="$(date +%s)"
deadline_epoch="$((start_epoch+1200))"
remote=/scratch/td2248/projects/single_market_verification_v1
source="$remote/source"
results="$remote/results"
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
baseline=/scratch/td2248/projects/grid_resolution_120x9_v1/results/full/control_160x15
[[ ! -e "$results" ]] || { echo 'Refusing existing results'; exit 2; }
(cd "$source" && sha256sum -c source.sha256)
[[ "$(sha256sum "$inputs/bundle.json" | cut -d' ' -f1)" == 427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7 ]]
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
mkdir -p "$results" "$remote/cache/preflight" "$remote/cache/production"
run_arm() {
  mode="$1"; cache="$2"; until="$3"
  export NUMBA_CACHE_DIR=/work/cache
  export PYTHONPATH=/work/source:/work/source/source
  remaining="$((until-$(date +%s)))"
  [[ "$remaining" -gt 0 ]]
  timeout --signal=TERM --kill-after=5s "${remaining}s" apptainer exec \
    --bind "$frozen:$repo:ro,$source:/work/source:ro,$inputs:/work/inputs:ro,$results:/work/results:rw,$remote/cache/$cache:/work/cache:rw" \
    --pwd "$repo" /share/apps/images/ubuntu-24.04.4.sif \
    /share/apps/anaconda3/2025.06/bin/python /work/source/run_arm.py "$mode" \
    --reference-root "$repo" --bundle /work/inputs --out "/work/results/$mode" --deadline-epoch "$deadline_epoch" \
    > "$results/$mode.log" 2>&1
}
# Authentication and exact mocked loop before any lifecycle work. Cold production cache.
run_arm smoke preflight "$((start_epoch+180))"
run_arm full production "$deadline_epoch"
# No wait or model work: if the baseline is not complete, preserve arm outputs for later comparison.
if [[ -e "$baseline/completed.json" ]]; then
  export PYTHONPATH=
  remaining="$((deadline_epoch-$(date +%s)))"
  [[ "$remaining" -gt 0 ]]
  timeout --signal=TERM --kill-after=5s "${remaining}s" /share/apps/anaconda3/2025.06/bin/python \
    "$source/compare.py" "$baseline" "$results/full" "$results/comparison.json"
else
  printf '{"status":"awaiting_baseline","baseline_job":18879780}\n' > "$results/comparison_pending.json"
fi
