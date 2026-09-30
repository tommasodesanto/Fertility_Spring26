#!/usr/bin/env bash
# Reviewed one-core validation batch: 2 component fixtures, then 2+2 fixed-price solves.
# No calibration, correction experiment or automatic retry. Failed step stops.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:40:00
#SBATCH --job-name=refactor_fp_pair
set -euo pipefail
LAB=/scratch/td2248/projects/publication_refactor_20260929
REF=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
ROOT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
SCALAR="$LAB/scalar_src_gridfix"
INDEXED="$LAB/indexed_src_gridfix"
export LAB BUNDLE_SHA=427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7
OUT="$LAB/results/fp_pair_gridfix_${SLURM_JOB_ID}"
mkdir "$OUT"
T0=$(date +%s)
CURRENT=preflight
finish() { rc=$?; printf '{"rc":%s,"last_phase":"%s","elapsed_seconds":%s,"max_lifecycle_solves":4}\n' "$rc" "$CURRENT" "$(( $(date +%s)-T0 ))" > "$OUT/completed.json"; }
trap finish EXIT
printf '{"maximum_lifecycle_solves":4,"components_cap_seconds":300,"each_fixed_price_phase_cap_seconds":900,"slurm_total_cap_seconds":2400,"estimated_seconds":650,"threads":1,"memory_gib":12,"stop":"first failure; no retries","economic_changes":"none; borrowing remains reference mode"}\n' > "$OUT/plan.json"
(cd "$SCALAR" && sha256sum -c --quiet SHA256SUMS)
(cd "$INDEXED" && sha256sum -c --quiet SHA256SUMS)
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
CURRENT=indexed_components
timeout -k 30 300 apptainer exec --bind "$REF:$ROOT:ro" --bind "$SCALAR:/lab:ro" --bind "$INDEXED:/idx:ro" --bind "$LAB:$LAB" --pwd /lab \
  --env PYTHONPATH=/lab:/idx,NUMBA_CACHE_DIR="$OUT/component_cache",REFACTOR_ROOT="$ROOT",REFACTOR_EXPORT="$LAB/export_v1",REFACTOR_BUNDLE_SHA="$BUNDLE_SHA",REFACTOR_INDEXED_STAGE=/idx/refactor_lab/engine,REFACTOR_INDEXED_MODULE=labidx.engine.kernels \
  /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python -m pytest /lab/refactor_lab/tests -k transformed_stage_kernel -q -p no:cacheprovider --junitxml="$OUT/indexed_components.xml" > "$OUT/indexed_components.log" 2>&1
CURRENT=scalar_fixed_price
LAB_SRC="$SCALAR" OUT_OVERRIDE="$OUT/scalar" PHASE=fixed-price MODE=torch bash "$SCALAR/refactor_lab/verify.sh" > "$OUT/scalar_driver.log" 2>&1
CURRENT=indexed_fixed_price
LAB_SRC="$INDEXED" OUT_OVERRIDE="$OUT/indexed" PHASE=fixed-price MODE=torch bash "$INDEXED/refactor_lab/verify.sh" > "$OUT/indexed_driver.log" 2>&1
CURRENT=all_passed
