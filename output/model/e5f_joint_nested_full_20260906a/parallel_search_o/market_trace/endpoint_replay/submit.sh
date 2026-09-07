#!/bin/bash
#SBATCH --job-name=joint_tax_endpoints
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/market_trace/endpoint_replay/slurm_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/market_trace/endpoint_replay/slurm_%j.err
set -euo pipefail
module purge
module load anaconda3/2025.06
export E5F_DIAGNOSTIC_SOURCE_ROOT=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o
export TMPDIR=/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/tmp
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/market_trace/endpoint_replay
python3 - <<'VERIFY'
import hashlib,pathlib
assert hashlib.sha256(pathlib.Path('diagnose_e5f_joint_policy_market.py').read_bytes()).hexdigest() == '67333dfb9a6b3a7dc988019b172fd07cd267a4e77d676e18c7945359846db769'
VERIFY
python3 -u diagnose_e5f_joint_policy_market.py --checkpoint /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/equilibrium_path/property-tax-2pct-no-rebate/date_2047/dated_state.pkl.gz --checkpoint-sha256 bd4c511443ea1cbc93133333c2b847c26cfec0b6eaf0756ef02079a3a95e1137 --progress /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/equilibrium_path/property-tax-2pct-no-rebate/policy_path_progress.csv --outdir /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/market_trace/endpoint_replay/results --endpoint-trace /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907o/output/model/joint_nested_overnight/market_trace/results/latest.json
