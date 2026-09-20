#!/bin/bash
#SBATCH --job-name=acs27-first-birth-housing
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/first_birth_housing_%j.out
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/first_birth_housing_%j.err
set -euo pipefail
ROOT="${KLEVEN_ROOT:-/scratch/td2248/projects/kleven_acs_pilot_20260917}"
SCRIPT="${ROOT}/code/empirical/acs/kleven_pseudo/run_first_birth_housing.R"
OUTDIR="${OUTDIR:-${ROOT}/output/kleven_acs_pilot/first_birth_housing_20260920}"
mkdir -p "${OUTDIR}"
exec > >(tee -a "${OUTDIR}/launcher.log") 2>&1
trap 's=$?; echo "launcher_exit=${s} line=${LINENO} time=$(date -Is)" >> "${OUTDIR}/launcher.log"; exit ${s}' EXIT
echo "startup time=$(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none}"
module load r/4.5.1
export R_LIBS_USER="${ROOT}/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
command -v Rscript
Rscript -e 'stopifnot(requireNamespace("data.table",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE),requireNamespace("fixest",quietly=TRUE)); cat("R_DEPENDENCY_SMOKE_PASS\\n")'
test -s "${SCRIPT}"
test -s "${ROOT}/overnight_ne_benchmark/ne_housing_v5/cps_acs_pseudo-panel.rds"
test -s "${ROOT}/output/kleven_acs_pilot/source_audit_extract27_20260919/ne_extract27_housing_key_packet.rds"
test -s "${ROOT}/output/kleven_acs_pilot/source_audit_extract27_20260919/source_audit_manifest.json"
export KLEVEN_ROOT="${ROOT}" ESTIMATOR_FILE="${ROOT}/code/empirical/acs/kleven_pseudo/estimate_first_birth_housing.R" OUTDIR
export RUN_ESTIMATION="${RUN_ESTIMATION:-1}"
export SECOND_BIRTH_DIAG_DIR="${ROOT}/output/kleven_acs_pilot/second_birth_proxy_diagnostic_20260920"
export CANDIDATE_OUTDIR="${ROOT}/output/kleven_acs_pilot/second_birth_candidate_support_20260920"
echo "preflight_pass run_estimation=${RUN_ESTIMATION} time=$(date -Is)"
OUTDIR="${CANDIDATE_OUTDIR}" Rscript "${ROOT}/code/empirical/acs/kleven_pseudo/second_birth_candidate_support.R"
Rscript "${SCRIPT}"
