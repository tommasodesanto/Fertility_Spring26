#!/bin/bash
# One approved allocation for the second-birth housing diagnostic.  It consumes
# the completed exact-support checkpoint and runs three declared specifications.
# No source reload or rematching occurs here.
#SBATCH --job-name=acs27-second-birth-housing
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/second_birth_housing_%j.out
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/second_birth_housing_%j.err

set -euo pipefail
ROOT="${KLEVEN_ROOT:-/scratch/td2248/projects/kleven_acs_pilot_20260917}"
CODE="${ROOT}/code/empirical/acs/kleven_pseudo"
MATCHED_DIR="${SECOND_BIRTH_MATCHED_DIR:-${ROOT}/output/kleven_acs_pilot/second_birth_transformed_support_18080547_corrected}"
OUTDIR="${SECOND_BIRTH_HOUSING_OUTDIR:-${ROOT}/output/kleven_acs_pilot/second_birth_housing_20260920}"
mkdir -p "${OUTDIR}"
exec > >(tee -a "${OUTDIR}/launcher.log") 2>&1
trap 's=$?; echo "launcher_exit=${s} line=${LINENO} time=$(date -Is)" >> "${OUTDIR}/launcher.log"; exit ${s}' EXIT

echo "startup time=$(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none}"
module load r/4.5.1
export R_LIBS_USER="${ROOT}/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript -e 'stopifnot(requireNamespace("data.table", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE), requireNamespace("fixest", quietly=TRUE)); cat("R_DEPENDENCY_SMOKE_PASS\n")'

# Verify the complete source closure before reading any checkpoint or fitting.
for f in \
  "${CODE}/run_second_birth_housing.R" \
  "${CODE}/second_birth_housing.R" \
  "${CODE}/build_second_birth_proxy.R" \
  "${CODE}/second_birth_matching.R" \
  "${MATCHED_DIR}/proxy_checkpoint.rds" \
  "${MATCHED_DIR}/matched_checkpoint.rds" \
  "${MATCHED_DIR}/packet_with_author_cells.rds"; do
  test -s "${f}"
done
echo "DEPENDENCY_CLOSURE_PASS"

cd "${CODE}"
KLEVEN_ROOT="${ROOT}" \
  SECOND_BIRTH_MATCHED_DIR="${MATCHED_DIR}" \
  SECOND_BIRTH_HOUSING_OUTDIR="${OUTDIR}" \
  Rscript run_second_birth_housing.R
