#!/bin/bash
#SBATCH --job-name=e5fnest_long
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=384G
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_nested_long_mpl_${SLURM_JOB_ID}"
cd "${SLURM_SUBMIT_DIR:?}"
: "${E5F_JOINT_MODE:?smoke or search required}"
: "${E5F_JOINT_CONTRACT:?required contract path}"
: "${E5F_JOINT_CONTRACT_SHA256:?required contract hash}"
CONTRACT="$E5F_JOINT_CONTRACT"
SHA="$E5F_JOINT_CONTRACT_SHA256"
printf '%s  %s\n' "$SHA" "$CONTRACT" | sha256sum --check -
python3 code/model/tools/test_e5f_exhaustive_saving.py
# A second process must also load the compiled kernel cache correctly.
python3 code/model/tools/test_e5f_exhaustive_saving.py
python3 code/model/tools/test_e5f_joint_nested_full.py
PYTHONPATH=code/model:code/model/tools python3 code/model/tools/test_e5f_joint_nested_integration.py
python3 code/model/tools/test_e5f_joint_nested_long_search.py
if [[ "$E5F_JOINT_MODE" == smoke ]]; then
 python3 -u code/model/tools/run_e5f_joint_nested_full_smoke.py \
  --checkpoint /scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/output/model/independent_numerical_smoke/dated_state.pkl \
  --output output/model/joint_nested_overnight/default_off_reference --reference
fi
if [[ "$E5F_JOINT_MODE" == policy-smoke ]]; then
 : "${E5F_JOINT_SELECTED_SUMMARY:?verified calibration summary required}"
 printf '%s  %s\n' "$SHA" "$CONTRACT" | sha256sum --check -
 exec python3 -u code/model/tools/run_e5f_joint_nested_finalize.py \
  --selected-summary "$E5F_JOINT_SELECTED_SUMMARY" \
  --outdir output/model/joint_nested_overnight/policy_loop_smoke \
  --contract "$CONTRACT" --smoke
fi
exec python3 -u code/model/tools/run_e5f_joint_nested_long_search.py --contract "$CONTRACT" --contract-sha256 "$SHA" --mode "$E5F_JOINT_MODE"
