#!/bin/bash
# One reviewed allocation for the short-window housing sensitivity and the
# transformed second-birth source-cell diagnostic.  The stages have separate
# output directories and are deliberately independent: a short-window support
# stop must not suppress the source-cell receipt.
#SBATCH --job-name=acs27-housing-second-support
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/housing_second_support_%j.out
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/housing_second_support_%j.err

set -euo pipefail
ROOT="${KLEVEN_ROOT:-/scratch/td2248/projects/kleven_acs_pilot_20260917}"
CODE="${ROOT}/code/empirical/acs/kleven_pseudo"
SECOND_OUT="${SECOND_BIRTH_OUTDIR:-${ROOT}/output/kleven_acs_pilot/second_birth_transformed_support_20260920}"
SHORT_OUT="${SHORT_WINDOW_OUTDIR:-${ROOT}/output/kleven_acs_pilot/first_birth_short_window_20260920}"
mkdir -p "${SECOND_OUT}" "${SHORT_OUT}"
exec > >(tee -a "${SHORT_OUT}/bundle_launcher.log") 2>&1
trap 's=$?; echo "bundle_exit=${s} line=${LINENO} time=$(date -Is)" >> "${SHORT_OUT}/bundle_launcher.log"; exit ${s}' EXIT

echo "startup time=$(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none}"
module load r/4.5.1
export R_LIBS_USER="${ROOT}/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript -e 'stopifnot(requireNamespace("data.table",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE),requireNamespace("fixest",quietly=TRUE)); cat("R_DEPENDENCY_SMOKE_PASS\n")'
test -s "${ROOT}/output/kleven_acs_pilot/source_audit_extract27_20260919/ne_extract27_housing_key_packet.rds"
test -s "${ROOT}/overnight_ne_benchmark/ne_acs.RData"
test -s "${CODE}/run_first_birth_short_window.R"
test -s "${CODE}/run_second_birth_transformed_support.R"
test -s "${CODE}/estimate_first_birth_housing.R"
test -s "${CODE}/build_second_birth_proxy.R"
test -s "${CODE}/second_birth_matching.R"
cd "${CODE}"

set +e
KLEVEN_ROOT="${ROOT}" SECOND_BIRTH_PACKET="${ROOT}/output/kleven_acs_pilot/source_audit_extract27_20260919/ne_extract27_housing_key_packet.rds" \
  SECOND_BIRTH_PREPARED="${ROOT}/overnight_ne_benchmark/ne_acs.RData" SECOND_BIRTH_OUTDIR="${SECOND_OUT}" \
  Rscript run_second_birth_transformed_support.R
second_status=$?
set -e
echo "second_birth_status=${second_status} time=$(date -Is)"

set +e
KLEVEN_ROOT="${ROOT}" OUTDIR="${SHORT_OUT}" ESTIMATOR_FILE="${CODE}/estimate_first_birth_housing.R" \
  Rscript run_first_birth_short_window.R
short_status=$?
set -e
echo "short_window_status=${short_status} time=$(date -Is)"

cat > "${SHORT_OUT}/bundle_stage_status.txt" <<EOF
job=${SLURM_JOB_ID:-none}
second_birth_status=${second_status}
short_window_status=${short_status}
second_birth_outdir=${SECOND_OUT}
short_window_outdir=${SHORT_OUT}
EOF

if [ "${second_status}" -ne 0 ] || [ "${short_status}" -ne 0 ]; then
  echo "bundle_status=FAILED_ONE_OR_MORE_STAGES"
  exit 1
fi
echo "bundle_status=COMPLETE"
