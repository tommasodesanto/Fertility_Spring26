#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:50:00
#SBATCH --job-name=refactor_ge_pair
set -euo pipefail
export LAB=/scratch/td2248/projects/publication_refactor_20260929
export LAB_SRC="$LAB/indexed_src_gridfix"
export BUNDLE_SHA=427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
BASE="$LAB/results/ge_pair_${SLURM_JOB_ID}"
mkdir "$BASE"
T0=$(date +%s); CURRENT=preflight
finish() { rc=$?; printf '{"rc":%s,"last_phase":"%s","elapsed_seconds":%s}\n' "$rc" "$CURRENT" "$(( $(date +%s)-T0 ))" > "$BASE/completed.json"; }
trap finish EXIT
printf '{"engines":"indexed lab and original frozen engine","start_factor":1.05,"threads":1,"memory_gib":12,"max_lifecycle_per_engine":18,"solve_seconds_per_engine":900,"phase_seconds":2700,"scheduler_seconds":3000,"estimated_total_seconds":850,"economic_changes":[],"renewal_threshold":"1e-6 difference from reference; diagnostic, not absolute renewal certification","stop":"first failure; no retries"}\n' > "$BASE/plan.json"
(cd "$LAB_SRC" && sha256sum -c --quiet SHA256SUMS)
echo "141deab226f35abb1981cf0473edf05709b67b768ed80eb0fec7d841366700e8  $LAB/strict_compare_v1.py" | sha256sum -c --quiet -
module load anaconda3/2025.06
CURRENT=matched_ge
OUT_OVERRIDE="$BASE/verify" PHASE=ge MODE=torch PRICE_FACTOR=1.05 RENEWAL_TOL=1e-6 bash "$LAB_SRC/refactor_lab/verify.sh" > "$BASE/driver.log" 2>&1
CURRENT=strict_all_array_comparison
timeout -k 30 120 /share/apps/anaconda3/2025.06/bin/python "$LAB/strict_compare_v1.py" --strict-paths --lab "$BASE/verify/lab/ge/solution_arrays.npz" --reference "$BASE/verify/old/solution_arrays.npz" --out "$BASE/strict_comparison.json" > "$BASE/strict_comparison.log" 2>&1
CURRENT=all_passed
