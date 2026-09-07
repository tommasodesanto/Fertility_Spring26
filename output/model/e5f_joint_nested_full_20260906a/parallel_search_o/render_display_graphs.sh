#!/bin/bash
#SBATCH --job-name=joint_review_graphs
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/display_graphs/slurm_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/display_graphs/slurm_%j.err
set -euo pipefail
module purge
module load anaconda3/2025.06
export E5F_DIAGNOSTIC_SOURCE_ROOT=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o
export TMPDIR=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/tmp
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/display_graphs
python3 - <<'VERIFY'
import hashlib,pathlib
assert hashlib.sha256(pathlib.Path('render_display_graphs.py').read_bytes()).hexdigest() == 'e674a4eff2c83ac1bcd330a2c3fc87d99bf5cf106d05af8c5f554e6d8510c969'
VERIFY
python3 render_display_graphs.py /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907m/output/model/joint_nested_overnight/smoke/smoke_histories/task_004 /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/display_graphs/results
