#!/bin/bash
#SBATCH --job-name=e5fnest_anchor
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_nested_anchor_mpl_${SLURM_JOB_ID}"
cd "${SLURM_SUBMIT_DIR:?}"
PLAN=output/model/joint_calibration_anchor/plan.json
SHA=8309bc530c648a11eb9ea4d11d435762b5caf906a340cb9ba054d21c07d3f3e2
python3 -u code/model/tools/run_e5f_joint_overnight_case.py --plan "$PLAN" --plan-sha256 "$SHA" --case-id 1 > output/model/joint_calibration_anchor/case_001.log 2>&1 &
PID1=$!
python3 -u code/model/tools/run_e5f_joint_overnight_case.py --plan "$PLAN" --plan-sha256 "$SHA" --case-id 2 > output/model/joint_calibration_anchor/case_002.log 2>&1 &
PID2=$!
RC1=0
RC2=0
wait "$PID1" || RC1=$?
wait "$PID2" || RC2=$?
if [[ "$RC1" != 0 || "$RC2" != 0 ]]; then exit 1; fi
