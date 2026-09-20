#!/usr/bin/env bash
# One native equilibrium/normalization pilot at fixed structural parameters.
# Stage the adapter, constructor, candidate and absolute-path plan before launch.
set -euo pipefail
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
ssh torch "sbatch --account=torch_pr_570_general --job-name=e5f_income_stationary --cpus-per-task=1 --mem=24G --time=01:00:00 --output=$experiment/income_stationary_%j.out" <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1 NUMBA_DISABLE_JIT=0
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
export NUMBA_CACHE_DIR="$experiment/numba_cache"
plan="$experiment/earnings_candidate/calibration_plan.remote.json"
output="$experiment/income_stationary_${SLURM_JOB_ID}"
cd "$experiment"
# Exact native loop: one fixed parameter point, at most eight fresh stationary
# solves to normalize completed fertility; 1800s native / 2100s scored wrapper.
# Prior native evaluation ~420s; changed income may need more normalization roots.
# Every numerical/source/target gate remains active; first failure stops the pilot.
python -B "$experiment/code/model/tools/run_e5f_income_candidate_calibration.py" --mode prepare-only --plan "$plan" --output "${output}_preflight"
set +e
timeout --signal=TERM --kill-after=30s 2400s python -B "$experiment/code/model/tools/run_e5f_income_candidate_calibration.py" --mode pilot --plan "$plan" --output "$output"
result=$?
printf '{"job_id":"%s","exit_code":%d,"scope":"fixed structural parameters; fresh income-candidate equilibrium and fit; no structural search"}\n' "$SLURM_JOB_ID" "$result" > "${output}_job_receipt.json"
exit "$result"
SBATCH
