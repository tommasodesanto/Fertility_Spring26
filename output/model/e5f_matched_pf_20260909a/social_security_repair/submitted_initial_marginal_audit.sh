#!/bin/bash
#SBATCH --job-name=e5f_initial_paygo
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_paygo_20260911a/audit_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_paygo_20260911a/audit_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=1 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_paygo_20260911a
export MPLCONFIGDIR="$PWD/mpl"
python - <<'CHECK'
import hashlib,json
from pathlib import Path
for n,pin in json.loads(Path('manifest.json').read_text())['files'].items():
    assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==pin,n
CHECK
python -m unittest test_e5f_stationary_paygo -v
python audit_e5f_stationary_paygo.py --contract /scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/contract.json --contract-sha256 13ac77b34c5b71094bf98461cb641260b7b32d130d1b88c599003132598d54a3 --source-root /scratch/td2248/projects/Fertility_Spring26_preference_shape_20260910b --output summary.json
