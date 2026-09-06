#!/bin/bash
#SBATCH --job-name=e5fnest_long
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=144G
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_nested_long_mpl_${SLURM_JOB_ID}"
cd "${SLURM_SUBMIT_DIR:?}"
: "${E5F_JOINT_MODE:?smoke or search required}"
CONTRACT=output/model/joint_nested_overnight/contract.json
SHA=416477db8ce66d22a6017aee24fd8a2a2d974c3fcf87bbed6bfe8f6f673c48ab
python3 code/model/tools/test_e5f_joint_nested_full.py
PYTHONPATH=code/model:code/model/tools python3 code/model/tools/test_e5f_joint_nested_integration.py
python3 code/model/tools/test_e5f_joint_nested_long_search.py
if [[ "$E5F_JOINT_MODE" == smoke ]]; then
 python3 -u code/model/tools/run_e5f_joint_nested_full_smoke.py \
  --checkpoint /scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/output/model/independent_numerical_smoke/dated_state.pkl \
  --output output/model/joint_nested_overnight/default_off_reference --reference
fi
exec python3 -u code/model/tools/run_e5f_joint_nested_long_search.py --contract "$CONTRACT" --contract-sha256 "$SHA" --mode "$E5F_JOINT_MODE"
