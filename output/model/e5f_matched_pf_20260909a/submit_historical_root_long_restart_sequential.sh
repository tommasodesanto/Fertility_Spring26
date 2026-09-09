#!/bin/bash
#SBATCH --job-name=e5f_long_hist_root_seq
#SBATCH --output=logs/historical_root_long_restart_sequential_%j.out
#SBATCH --error=logs/historical_root_long_restart_sequential_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Collector 17279797 completed and its receipt was verified before submission.
# Fill the contract's jacobian_packet_sha256 from the verified long receipt,
# then replace the contract hash below with SHA256 of that finalized contract.
# Four complete 28-date paths = 224 Bellman solves, including fresh final replay.
# Latest mean 783 seconds/path: about 52 minutes; internal watchdog 70 minutes.
set -euo pipefail
CONTRACT_SHA256='bd2b699be6d1d28d347cd6ab67481f59070a955e5cd839edc029500b4fb4a515'
if [[ ! "$CONTRACT_SHA256" =~ ^[0-9a-f]{64}$ ]]; then
    echo 'BLOCKED: pending verified long Jacobian and finalized contract hash.' >&2
    exit 2
fi
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/mpl_root_long_restart_sequential
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/numba_root_long_restart_sequential
export TMPDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/tmp_root_long_restart_sequential
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$TMPDIR"
printf '%s  %s\n' "$CONTRACT_SHA256" historical_root_long_restart_sequential_contract.json | sha256sum --check --status
python code/model/tools/run_e5f_matched_pf_historical_root.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/historical_root_long_restart_sequential_contract.json --contract-sha256 "$CONTRACT_SHA256" --arm sequential --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/historical_root_long_restart_01/sequential
