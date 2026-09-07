#!/bin/bash
#SBATCH --job-name=e5f_two_shock
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=1-2
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="code/model:code/model/tools"
cd "${SLURM_SUBMIT_DIR:?}"
: "${E5F_TWO_SHOCK_PLAN:?}"
: "${E5F_TWO_SHOCK_PLAN_SHA256:?}"
# The two identical cases are the full historical-loop smoke, not a search.
python3 code/model/tools/test_e5f_two_shock_calibration.py
python3 code/model/tools/test_e5f_two_shock_choice.py
python3 code/model/tools/test_e5f_two_shock_choice.py
exec python3 -u code/model/tools/run_e5f_bounded_calibration_refinement.py \
 --plan "$E5F_TWO_SHOCK_PLAN" --plan-sha256 "$E5F_TWO_SHOCK_PLAN_SHA256" \
 --case-id "${SLURM_ARRAY_TASK_ID:?}"
