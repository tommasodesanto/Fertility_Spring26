#!/bin/bash
#SBATCH --job-name=e5f_active_optimizer
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_active_optimizer_20260911a/diagnostic_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_active_optimizer_20260911a/diagnostic_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_paygo_parenthood_1a5e88b7
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
python - <<'PY'
import hashlib,json
from pathlib import Path
for name,pin in json.loads(Path('source_pins.json').read_text()).items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==pin,name
p=Path('/scratch/td2248/projects/Fertility_Spring26_active_optimizer_20260911a/active_optimizer_diagnostic.py')
assert hashlib.sha256(p.read_bytes()).hexdigest()=='a4e975c81679d11578c52627c0a9219c34fcc80dc08336e103e0e23719c5e552'
PY
python /scratch/td2248/projects/Fertility_Spring26_active_optimizer_20260911a/active_optimizer_diagnostic.py --repo-root "$PWD" --output-dir /scratch/td2248/projects/Fertility_Spring26_active_optimizer_20260911a/output
