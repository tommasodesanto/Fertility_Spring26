#!/usr/bin/env bash
#SBATCH --job-name=surv_exit_exp
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
base=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1
control=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/estate_receiver_probe_20260925_v1/results/v2/run/control
out=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/survival_rationale_review_20260925/household_mapping
export PYTHONPATH="$base/source/code/model/tools:$base/source/code/model${PYTHONPATH:+:$PYTHONPATH}"
python3 "$out/measure_checkpoint_exposure.py" \
  --packet "$control/initial_state.pkl.gz" \
  --receipt "$control/receipt.json" \
  --lifecycle "$control/lifecycle_2023.csv" \
  --schedule "$out/candidate_schedule.csv" \
  --output "$out/checkpoint_exposure.json"
