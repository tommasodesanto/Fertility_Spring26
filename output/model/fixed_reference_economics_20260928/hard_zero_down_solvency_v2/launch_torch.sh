#!/usr/bin/env bash
#SBATCH --job-name=hard_solvency_pe
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/hard_zero_down_solvency_v2/logs/%x-%j.out
set -euo pipefail
remote=/scratch/td2248/projects/hard_zero_down_solvency_v2
strict_remote=/scratch/td2248/projects/strict_purchase_sandbox_v1
v2=/scratch/td2248/projects/normalized_calibration_v2
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/hard_zero_down_solvency_v2
v2packet=output/model/fixed_reference_economics_20260928/normalized_calibration_v2
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
out="$remote/results/onecase"
mkdir "$out" || { echo 'Refusing existing results'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal() {
  code=$?
  trap - EXIT
  "$python" - "$out" "$code" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1], 'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$code"
}
trap terminal EXIT
trap 'exit 143' TERM
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$v2/source/$v2packet:$repo/$v2packet:ro")
while IFS= read -r rel; do binds+=(--bind "$v2/source/$rel:$repo/$rel:ro"); done < <("$python" - "$v2/inventory.json" "$v2packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if not rel.startswith(sys.argv[2]+'/'): print(rel)
PY
)
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
binds+=(--bind "$strict_remote/source/code/model/experiments/strict_purchase_sandbox:$repo/code/model/experiments/strict_purchase_sandbox:ro")
binds+=(--bind "$remote/solvency_source:$repo/code/model/experiments/strict_purchase_solvency/source:ro")
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
timeout --signal=TERM --kill-after=10s 1180s apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/run.py" --out /work/results/run > "$out/run.log" 2>&1
