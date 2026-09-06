#!/bin/bash
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_reporting_${SLURM_JOB_ID}"
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906e
CHECKPOINT=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906c/output/model/joint_nested_overnight/smoke/smoke_anchor/task_001/dated_state.pkl.gz
SHA=35b7624f605b2ed75f47c025a751805032901c953efe72637ebb9a440c130ac7
python3 -u code/model/tools/run_e5f_joint_nested_reporting_check.py --model-root /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906d --bundle-sha256 2c4adcdec5b85cff39b4c5d1466e6224db1dddf13052fe66c8a6be86e9a32951 --checkpoint "$CHECKPOINT" --checkpoint-sha256 "$SHA" --outdir output/model/joint_nested_overnight/reporting_before
python3 -u code/model/tools/run_e5f_joint_nested_reporting_check.py --model-root /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906e --bundle-sha256 85450db0d7611f7206fba933a74c0f962c18990917a926f9bb2888057494ff39 --checkpoint "$CHECKPOINT" --checkpoint-sha256 "$SHA" --outdir output/model/joint_nested_overnight/reporting_after --reference output/model/joint_nested_overnight/reporting_before/reporting_check.json
