#!/usr/bin/env bash
# Upload the bounded credit diagnostic beside immutable v5 and submit it.
#SBATCH --job-name=cescredit
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v5/logs/cescredit-%j.out
set -euo pipefail

remote=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v5
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
sha=26bd0fbb6901f9080d4e15472029c48847f117b21cb95f6d0de0be4cc3976e85

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  module load anaconda3/2025.06
  export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
  out="$remote/results/credit_diagnostic_${SLURM_JOB_ID}"
  mkdir "$out" "$out/numba_cache" "$out/matplotlib"
  write_terminal() {
    local s=$?
    "$python" - "$out" "$s" <<'PY'
import json,sys,time
from pathlib import Path
Path(sys.argv[1],"launcher_terminal.json").write_text(json.dumps(dict(exit_code=int(sys.argv[2]),finished_epoch=time.time()),indent=2)+"\n")
PY
    exit "$s"
  }
  trap write_terminal EXIT
  export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
  container_env=(--env CES_NORMALIZED_SHARES_STAGED_CONTEXT=1 --env "CES_CREDIT_CHAIN_SHA256=$sha" --env "PYTHONPATH=$repo/code/model")
  binds=(--bind "$remote:/work/deployment:ro" --bind "$remote/source:$repo:ro" --bind "$out:/work/results:rw")
  apptainer exec "${container_env[@]}" "${binds[@]}" --pwd "$repo" \
    /share/apps/images/ubuntu-24.04.4.sif "$python" /work/deployment/verify_stage.py --container
  apptainer exec "${container_env[@]}" "${binds[@]}" --pwd "$repo" \
    /share/apps/images/ubuntu-24.04.4.sif "$python" /work/deployment/followup_tools/credit_diagnostic.py \
    --prepare-only --out /work/results/preflight
  timeout --signal=TERM --kill-after=10s 1700s apptainer exec "${container_env[@]}" "${binds[@]}" --pwd "$repo" \
    /share/apps/images/ubuntu-24.04.4.sif "$python" /work/deployment/followup_tools/credit_diagnostic.py \
    --out /work/results/credit_diagnostic
  exit 0
fi

here=$(cd "$(dirname "$0")" && pwd)
chain="$repo/output/model/experiments/ces_normalized_shares/overnight_v1/final_results/chain_1.json"
actual=$(shasum -a 256 "$chain" | cut -d ' ' -f 1)
[[ "$actual" == "$sha" ]] || { echo "chain-1 input SHA-256 mismatch" >&2; exit 2; }
ssh -o BatchMode=yes -o ConnectTimeout=12 torch "mkdir -p '$remote/followup_tools' '$remote/results'"
rsync --checksum -e 'ssh -o BatchMode=yes' "$here/credit_diagnostic.py" \
  "torch:$remote/followup_tools/credit_diagnostic.py"
rsync --checksum -e 'ssh -o BatchMode=yes' "$chain" \
  "torch:$remote/followup_tools/credit_reference_chain1.json"
rsync --checksum -e 'ssh -o BatchMode=yes' "$here/credit_diagnostic.sh" \
  "torch:$remote/followup_tools/credit_diagnostic.sh"
ssh -o BatchMode=yes torch "'$python' '$remote/verify_stage.py' --host && sbatch '$remote/followup_tools/credit_diagnostic.sh'"
