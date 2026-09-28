#!/usr/bin/env bash
#SBATCH --job-name=e5f_price_diagnostic
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:45:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
diag="$stage/price_diagnostic_v1"
run="$stage/project/output/model/overnight_calibration_20260928/gated_v1"
hard_end=1790595900 #07:45 EDT, reserve morning review time
# Only after search dispatch has ended; final repeats use at most eight workers.
while (( $(date +%s) < 1790593200 )); do sleep 15; done
[ -f "$run/lead_search_approval.json" ] || exit 3
while :; do
  now=$(date +%s)
  (( now < hard_end - 600 )) || { echo 'Too late for bounded diagnosis'; exit 3; }
  active=$(timeout 10 srun --jobid=18687184 --overlap -N1 -n1 -c1 ps -u "$(id -u)" -o args= 2>/dev/null | awk '/--stage evaluate/ && /overnight_calibration_20260928\/gated_v1\/search\// && !/awk/ {n++} END {print n+0}') || active=unknown
  # If main allocation ended, Slurm reports no running job; no model workers remain there.
  if [ -z "$(squeue -j18687184 -h -o '%T')" ]; then active=0; fi
  if [[ "$active" =~ ^[0-9]+$ ]] && (( active <= 8 )); then break; fi
  sleep 15
done
printf '{"epoch":%s,"main_model_workers":%s,"diagnostic_workers":1,"maximum_diagnostic_GE_calls":2}\n' "$(date +%s)" "$active" > "$diag/concurrency_at_launch.json"
cd "$diag"
sha256sum -c approved_manifest.sha256
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$original/code/model/tools"
export NUMBA_CACHE_DIR="$stage/numba_cache"
remaining=$((hard_end - $(date +%s) - 10))
exec timeout --signal=TERM --kill-after=10 "$remaining" apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$diag/diagnose_e5f_evening_housing_failures.py" --plan "$diag/approved_plan.json" --plan-sha256 "$(cat "$diag/approved_plan.sha256")" --output "$diag/run_v1"
