#!/usr/bin/env bash
# Submit only after staging helpers and candidate.json to the independent Torch folder.
set -euo pipefail
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
ssh torch "sbatch --account=torch_pr_570_general --job-name=e5f_income_cohort --cpus-per-task=1 --mem=24G --time=01:00:00 --output=$experiment/income_%j.out" <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
native=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches
out="$experiment/income_cohort_${SLURM_JOB_ID}"
mkdir -p "$out"
export NUMBA_CACHE_DIR="$experiment/numba_cache"
args=(--checkpoint "$native/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz" --replay "$experiment/replay" --source-root "$native/final_night_20260913/corrected_initial_source_v2/code/model" --prior-input "$experiment/run" --candidate-json "$experiment/earnings_candidate/candidate.json" --output "$out")
# Estimate: two ~70s Bellman solves, 68 cohort age maps, standard diagnostics;
# allow 55 minutes, fail on the first numerical/contract violation. No refit.
python -B "$experiment/code/model/tools/run_e5f_native_income_cohort_diagnostic.py" --mode inspect "${args[@]}"
set +e
timeout --signal=TERM --kill-after=30s 3300s python -B "$experiment/code/model/tools/run_e5f_native_income_cohort_diagnostic.py" --mode candidate "${args[@]}"
result=$?
printf '{"job_id":"%s","exit_code":%d,"scope":"earnings cohort diagnostic; no recalibration"}\n' "$SLURM_JOB_ID" "$result" > "$out/job_receipt.json"
exit "$result"
SBATCH
