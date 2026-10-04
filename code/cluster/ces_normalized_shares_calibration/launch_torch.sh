#!/usr/bin/env bash
# Immutable CES stage. preflight is zero-solve; smoke/production require Slurm.
#SBATCH --job-name=cesshares
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3/logs/%x-%A_%a.out
set -euo pipefail
remote=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
mode=${1:-${CES_RUN_MODE:-production}}
[[ "$mode" =~ ^(preflight|smoke|production)$ ]] || { echo "invalid mode"; exit 2; }
if [[ "$mode" == preflight ]]; then chain=0; else chain=${SLURM_ARRAY_TASK_ID:?array 0..3 required}; fi
[[ "$chain" =~ ^[0-3]$ ]] || { echo "invalid chain"; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
"$python" "$remote/verify_stage.py" --host
starts_rel=output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json
starts_sha=$("$python" - "$remote/inventory.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['start_plan_sha256'])
PY
)
# Staged source is mounted last and read-only, including parent compatibility/observer files.
binds=(--bind "$remote:/work/deployment:ro" --bind "$remote/source:$repo:ro")
container_env=(--env CES_NORMALIZED_SHARES_STAGED_CONTEXT=1 --env "PYTHONPATH=$repo/code/model")
if [[ "$mode" == preflight ]]; then
  out="$remote/preflight/chain_0"; mkdir -p "$out"
  apptainer exec "${container_env[@]}" "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/verify_stage.py --container
  apptainer exec "${container_env[@]}" "${binds[@]}" --bind "$out:/work/results:rw" --pwd "$repo" "$image" "$python" "$repo/code/model/experiments/ces_normalized_shares/calibrate.py" --chain 0 --out /work/results/run --deadline-epoch "$(( $(date +%s)+3600 ))" --starts-file "$repo/$starts_rel" --starts-file-sha256 "$starts_sha" --preflight-evaluator > "$out/driver.log" 2>&1
  "$python" - "$out/run/completed.json" <<'PY'
import json,sys
assert json.load(open(sys.argv[1]))['status']=='evaluator_initialized_zero_solves'
PY
  exit 0
fi
start=$(date +%s); wall=21600; [[ "$mode" == smoke ]] && wall=5400; deadline=$((start+wall)); out="$remote/results/${mode}_chain_${chain}"
mkdir "$out" || { echo "refusing existing result directory"; exit 2; }; mkdir "$out/numba_cache" "$out/matplotlib"
"$python" - "$out" "$mode" "$chain" "$start" "$deadline" "$wall" "$remote" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
out,mode,chain,start,deadline,wall,remote=sys.argv[1:]
Path(out,'launcher_start.json').write_text(json.dumps(dict(mode=mode,chain=int(chain),start_epoch=int(start),deadline_epoch=int(deadline),wall_seconds=int(wall),cpus=1,memory_GiB=24,maximum_objective_calls=500,final_native_reserve_seconds=1800,stage_inventory_sha256=hashlib.sha256((Path(remote)/'inventory.json').read_bytes()).hexdigest(),no_auto_retry=True),indent=2)+'\n')
PY
trap 's=$?; "$python" - "$out" "$s" <<'"'"'PY'"'"'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],"launcher_terminal.json").write_text(json.dumps(dict(exit_code=int(sys.argv[2]),finished_epoch=time.time(),no_auto_retry=True),indent=2)+"\\n")
PY
exit "$s"' EXIT
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
opts=(--chain "$chain" --out /work/results/run --deadline-epoch "$deadline" --starts-file "$repo/$starts_rel" --starts-file-sha256 "$starts_sha")
[[ "$mode" == smoke ]] && opts+=(--smoke)
timeout --signal=TERM --kill-after=10s "$((deadline-$(date +%s)-15))s" apptainer exec "${container_env[@]}" "${binds[@]}" --bind "$out:/work/results:rw" --pwd "$repo" "$image" "$python" "$repo/code/model/experiments/ces_normalized_shares/calibrate.py" "${opts[@]}" > "$out/driver.log" 2>&1
