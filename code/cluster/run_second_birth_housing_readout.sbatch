#!/bin/bash
# Bounded readout allocation for the completed second-birth housing run.
# It reuses saved checkpoints, performs no rematching, and gates the common-
# cohort sensitivity before any diagnostic fit.
#SBATCH --job-name=acs27-second-birth-readout
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/second_birth_readout_%j.out
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/second_birth_readout_%j.err

set -euo pipefail
ROOT="${KLEVEN_ROOT:-/scratch/td2248/projects/kleven_acs_pilot_20260917}"
CODE="${ROOT}/code/empirical/acs/kleven_pseudo"
MATCHED_DIR="${SECOND_BIRTH_MATCHED_DIR:-${ROOT}/output/kleven_acs_pilot/second_birth_transformed_support_18080547_corrected}"
HOUSING_DIR="${SECOND_BIRTH_HOUSING_OUTDIR:-${ROOT}/output/kleven_acs_pilot/second_birth_housing_20260920_job_a218afc6}"
OUTDIR="${SECOND_BIRTH_READOUT_OUTDIR:-${ROOT}/output/kleven_acs_pilot/second_birth_housing_readout_20260920}"
mkdir -p "${OUTDIR}"
exec > >(tee -a "${OUTDIR}/launcher.log") 2>&1
trap 's=$?; echo "launcher_exit=${s} line=${LINENO} time=$(date -Is)" >> "${OUTDIR}/launcher.log"; exit ${s}' EXIT

echo "startup time=$(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none}"
module load r/4.5.1
export R_LIBS_USER="${ROOT}/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript -e 'stopifnot(requireNamespace("data.table",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE),requireNamespace("fixest",quietly=TRUE)); cat("R_DEPENDENCY_SMOKE_PASS\n")'
for f in \
  "${CODE}/run_second_birth_housing_readout.R" \
  "${CODE}/second_birth_housing.R" \
  "${MATCHED_DIR}/proxy_checkpoint.rds" \
  "${MATCHED_DIR}/matched_checkpoint.rds" \
  "${MATCHED_DIR}/packet_with_author_cells.rds" \
  "${HOUSING_DIR}/primary/support.csv" \
  "${HOUSING_DIR}/primary/curves.csv"; do test -s "${f}"; done
echo "READOUT_DEPENDENCY_CLOSURE_PASS"
cd "${CODE}"
KLEVEN_ROOT="${ROOT}" SECOND_BIRTH_MATCHED_DIR="${MATCHED_DIR}" \
  SECOND_BIRTH_HOUSING_OUTDIR="${HOUSING_DIR}" SECOND_BIRTH_READOUT_OUTDIR="${OUTDIR}" \
  Rscript run_second_birth_housing_readout.R
