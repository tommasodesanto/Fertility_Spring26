#!/usr/bin/env bash
#SBATCH --job-name=cesoldprice
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v5/logs/cesoldprice-%j.out
set -euo pipefail
remote=/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v5
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
out="$remote/results/historical_start_upperprice"
mkdir "$out"
mkdir "$out/numba_cache" "$out/matplotlib"
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
"$python" "$remote/verify_stage.py" --host
timeout --signal=TERM --kill-after=10s 1600s apptainer exec --env CES_NORMALIZED_SHARES_STAGED_CONTEXT=1 --env "PYTHONPATH=$repo/code/model" --bind "$remote:/work/deployment:ro" --bind "$remote/source:$repo:ro" --bind "$out:/work/results:rw" --pwd "$repo" /share/apps/images/ubuntu-24.04.4.sif "$python" /work/deployment/followup_tools/check_historical_upperprice.py > "$out/driver.log" 2>&1
