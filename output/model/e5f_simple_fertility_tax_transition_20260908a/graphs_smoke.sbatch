#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/graph-%A_%a.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/mpl
export TMPDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/tmp-${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$MPLCONFIGDIR"
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
python -c "import hashlib; from pathlib import Path; assert hashlib.sha256(Path('code/model/tools/run_e5f_simple_fertility_tax_graphs.py').read_bytes()).hexdigest() == 'f6d41ae3953daef48d146fa955f481f939261416a5567b06c1ad8bbf5f87fb31'"
case=tax1-equal-rebate
year=2023
python code/model/tools/run_e5f_simple_fertility_tax_graphs.py --transition-contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/contract.json --transition-contract-sha256 ce80b6ad241bec9556d2f3d3cdccb9f89e0ccfdd5b2e80ddcd167ed5156b5b68 --stage smoke --case "$case" --year "$year" --output-root /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/graphs/smoke
