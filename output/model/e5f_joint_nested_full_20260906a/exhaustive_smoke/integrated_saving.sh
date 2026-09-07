#!/bin/bash
#SBATCH --job-name=e5fnest_saving
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:25:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_nested_saving_${SLURM_JOB_ID}"
cd "${SLURM_SUBMIT_DIR:?}"
python3 code/model/tools/test_e5f_exhaustive_saving.py
python3 -u code/model/tools/run_e5f_joint_nested_full_smoke.py \
 --checkpoint /scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/output/model/independent_numerical_smoke/dated_state.pkl \
 --output output/model/saving_integration/default_off_reference --reference
python3 -u code/model/tools/run_e5f_joint_nested_reporting_check.py \
 --model-root /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907b \
 --bundle-sha256 733ccb1b975e55d5baf1a46d7733affe0b0272dd9e5cfc4d4e38e25a32d1e387 \
 --checkpoint /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906e/output/model/joint_nested_overnight/policy_loop_smoke/supply-plus-20/date_2027/dated_state.pkl.gz \
 --checkpoint-sha256 bae769e9f08434c3d21b84d7f4e01573e6fa79ddd8c4f3cf4fa65d77a570513e \
 --outdir output/model/saving_integration/comparison --saving-integration
