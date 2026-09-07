#!/bin/bash
#SBATCH --job-name=e5f_renter_verify
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=96G
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907j
export TMPDIR="$PWD/tmp" NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$TMPDIR/mpl"
python3 code/model/tools/test_e5f_exhaustive_saving.py
python3 code/model/tools/test_e5f_exhaustive_saving.py
python3 code/model/tools/test_e5f_joint_nested_long_search.py
python3 code/model/tools/test_e5f_joint_nested_full.py
PYTHONPATH=code/model:code/model/tools python3 code/model/tools/test_e5f_joint_nested_integration.py
python3 code/model/tools/test_e5f_small_mass_transport.py
python3 -u code/model/tools/run_e5f_joint_nested_full_smoke.py --checkpoint /scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905/output/model/independent_numerical_smoke/dated_state.pkl --output output/model/joint_nested_overnight/default_off_reference --reference
python3 -u verify_and_replay_case17.py > logs/case17_replay.log 2>&1 &
replay_pid=$!
trap 'kill "$replay_pid" 2>/dev/null || true' EXIT
E5F_JOINT_MODE=smoke E5F_JOINT_CONTRACT="$PWD/output/model/joint_nested_overnight/contract.json" E5F_JOINT_CONTRACT_SHA256="$(sha256sum output/model/joint_nested_overnight/contract.json | cut -d ' ' -f1)" bash code/cluster/submit_e5f_joint_nested_long.sh
wait "$replay_pid"
trap - EXIT
