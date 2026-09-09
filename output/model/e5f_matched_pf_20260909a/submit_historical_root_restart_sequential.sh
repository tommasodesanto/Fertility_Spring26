#!/bin/bash
#SBATCH --job-name=e5f_hist_restart_seq
#SBATCH --output=logs/historical_root_restart_sequential_%j.out
#SBATCH --error=logs/historical_root_restart_sequential_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:35:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# PREPARED ONLY: verify completion of 17279004 and all 12 target-fit replay rows.
# Fill all three restart receipt hashes in the contract, then replace the hash below
# with SHA256 of that finalized contract. Original short Jacobian remains pinned.
# Four complete 12-date paths = 96 Bellman solves, including fresh final replay.
# About 26 minutes at measured timing; internal watchdog 30 minutes.
set -euo pipefail
CONTRACT_SHA256='899c219bdf648b749be292a8ea66436d1b03d24f21e144b0dc17013b085d763c'
if [[ ! "$CONTRACT_SHA256" =~ ^[0-9a-f]{64}$ ]]; then
    echo 'BLOCKED: pending completed root receipts, verified replay and finalized contract hash.' >&2
    exit 2
fi
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/mpl_root_restart_sequential
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/numba_root_restart_sequential
export TMPDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/tmp_root_restart_sequential
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$TMPDIR"
printf '%s  %s\n' "$CONTRACT_SHA256" historical_root_restart_sequential_contract.json | sha256sum --check --status
python code/model/tools/run_e5f_matched_pf_historical_root.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/historical_root_restart_sequential_contract.json --contract-sha256 "$CONTRACT_SHA256" --arm sequential --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/historical_root_restart_01/sequential
