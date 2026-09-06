#!/bin/bash
#SBATCH --job-name=e5fnest_smoke
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:25:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_nested_mpl_${SLURM_JOB_ID}"
cd "${SLURM_SUBMIT_DIR:?}"
python3 code/model/tools/test_e5f_joint_nested_full.py
PYTHONPATH=code/model:code/model/tools python3 code/model/tools/test_e5f_joint_nested_integration.py
python3 -u code/model/tools/run_e5f_joint_nested_full_smoke.py \
 --checkpoint /scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/output/model/independent_numerical_smoke/dated_state.pkl \
 --output output/model/core_smoke --scale 2.0 --lambda 0.8
