#!/usr/bin/env bash
# NE-only roster/proxy support diagnostic. This stops before matching or housing estimation.
#SBATCH --job-name=acs27-proxy-diagnostic
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/second_birth_proxy_diagnostic_20260920/slurm_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/second_birth_proxy_diagnostic_20260920/slurm_%j.log
set -euo pipefail

ROOT=${KLEVEN_ROOT:-/scratch/td2248/projects/kleven_acs_pilot_20260917}
AUD=$ROOT/forensics/source_audit_extract27_20260919
OUT=$ROOT/output/kleven_acs_pilot/second_birth_proxy_diagnostic_20260920
mkdir -p "$OUT"
step=before_module
diag="$OUT/launcher_${SLURM_JOB_ID}.log"
exec > >(tee -a "$diag") 2>&1
trap 'rc=$?; printf "exit_step=%s rc=%s\n" "$step" "$rc" >&2' EXIT
printf 'job_id=%s host=%s cwd=%s\n' "$SLURM_JOB_ID" "$(hostname)" "$PWD"

step=module_load
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT"

step=runtime_smoke
command -v Rscript
Rscript --version
Rscript --vanilla -e 'pkgs <- c("data.table","jsonlite"); stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE))); cat("R_RUNTIME_SMOKE_OK\n")'

step=proxy_diagnostic
Rscript --vanilla "$AUD/run_second_birth_proxy_diagnostic.R"
step=complete
