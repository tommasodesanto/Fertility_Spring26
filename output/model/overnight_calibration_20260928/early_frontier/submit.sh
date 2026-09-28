#!/usr/bin/env bash
#SBATCH --job-name=e5f_early_frontier
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=00:45:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
diag="$stage/early_frontier_v1"
main="$stage/project/output/model/overnight_calibration_20260928/gated_v1/search"
cd "$diag"
sha256sum -c approved_manifest.sha256
module load anaconda3/2025.06
while (( $(date +%s) < 1790593200 )); do sleep 10; done
proof_args=()
while :; do
  (( $(date +%s) < 1790593800 )) || { echo 'Insufficient remaining diagnostic window'; exit 3; }
  state=$(squeue -j18687184 -h -o '%T' 2>/dev/null || true)
  if [ -n "$state" ]; then
    active=$(timeout 10 srun --jobid=18687184 --overlap -N1 -n1 -c1 ps -u "$(id -u)" -o args= 2>/dev/null | awk '/--stage evaluate/ && /overnight_calibration_20260928\/gated_v1\/search\// && !/awk/ {n++} END {print n+0}') || active=unknown
    if [[ "$active" =~ ^[0-9]+$ ]] && (( active <= 8 )); then break; fi
  else
    status=$(sacct -j18687184 -X -n -P --format=State,ExitCode | head -1)
    [ "$status" = 'COMPLETED|0:0' ] || { echo "Main job not certified complete: $status"; exit 4; }
    python - "$main/complete.json" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
p=Path(sys.argv[1]);x=json.loads(p.read_text());assert x['status']=='bounded_search_complete'
proof={'job_id':'18687184','state':'COMPLETED','exit_code':'0:0','verified_epoch':time.time(),'complete_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
Path('main_completion_proof.json').write_text(json.dumps(proof)+'\n')
PY
    proof_sha=$(sha256sum main_completion_proof.json | cut -d' ' -f1)
    proof_args=(--completion-proof "$diag/main_completion_proof.json" --completion-proof-sha256 "$proof_sha")
    active=0
    break
  fi
  sleep 10
done
printf '{"epoch":%s,"main_workers":%s,"frontier_workers_max":10,"other_diagnostic_workers_max":1,"combined_max":19}\n' "$(date +%s)" "$active" > concurrency_at_launch.json
unset NUMBA_DISABLE_JIT APPTAINERENV_NUMBA_DISABLE_JIT SINGULARITYENV_NUMBA_DISABLE_JIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$original/code/model/tools"
export NUMBA_CACHE_DIR="$stage/numba_cache"
remaining=$((1790595900 - $(date +%s) - 10))
exec timeout --signal=TERM --kill-after=10 "$remaining" apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$diag/run_e5f_early_fertility_frontier.py" --plan "$diag/approved_plan.json" --plan-sha256 "$(cat approved_plan.sha256)" --output "$diag/run_v1" "${proof_args[@]}"
