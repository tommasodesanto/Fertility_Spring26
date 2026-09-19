#!/usr/bin/env bash
# Submit from the staged experiment directory on Torch.
# Five full-household solves: two controls, then three financing treatments.
# Expected 15-25 minutes from prior fixed-price timings; hard budget 45 minutes.
# No calibration, price root, preference normalization, or population transition.
#SBATCH --account=torch_pr_570_general
#SBATCH --job-name=native_finance
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:55:00
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
native=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches
export NUMBA_CACHE_DIR="$experiment/numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"
python -B "$experiment/code/model/tools/run_e5f_native_financing_diagnostic.py" \
  --mode sequence \
  --checkpoint "$native/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz" \
  --replay "$experiment/replay" \
  --source-root "$native/final_night_20260913/corrected_initial_source_v2/code/model" \
  --output "$experiment/run" \
  --case-budget-seconds 600 --total-budget-seconds 2700
