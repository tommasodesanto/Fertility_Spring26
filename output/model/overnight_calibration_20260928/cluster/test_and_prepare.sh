#!/usr/bin/env bash
#SBATCH --job-name=e5f_night_preflight
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
sha256sum -c "$stage/refresh_manifest.sha256"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export REVIEWED_RECOVERY_SEARCH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4/run_e5f_utility_comparison_search.py"
# Synthetic controller fixtures only. No household solve is called.
export NUMBA_DISABLE_JIT=1
apptainer exec --bind "$stage/project:$original" --pwd "$original/code/model/tools" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python -m unittest -v test_e5f_night_calibration
unset NUMBA_DISABLE_JIT
apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python code/model/tools/prepare_e5f_night_contract.py \
 --source-root "$original" \
 --evening-contract "$original/output/model/evening_calibration_20260927/contract_v4/contract.json" \
 --anchor-case "$original/output/model/evening_calibration_20260927/gated_v2/search/initial_0347_block/case" \
 --anchor-repeat-one "$original/output/model/evening_calibration_20260927/gated_v2/search/repeat_0364_block/case" \
 --anchor-repeat-two "$original/output/model/evening_calibration_20260927/gated_v2/search/repeat_0365_block/case" \
 --start-epoch 1790563080 \
 --output "$original/output/model/overnight_calibration_20260928/contract_v1"
contract="$stage/project/output/model/overnight_calibration_20260928/contract_v1/contract.json"
export EXPECTED_E5F_NIGHT_SHA256=$(sha256sum "$contract" | awk '{print $1}')
apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python code/model/tools/run_e5f_night_calibration.py \
 --stage prepare --contract "$original/output/model/overnight_calibration_20260928/contract_v1/contract.json" \
 --output "$original/output/model/overnight_calibration_20260928/contract_v1/controller_preparation"
