#!/usr/bin/env bash
# Source audit for the verified IPUMS extract27 housing payload.
# The source digest was independently verified before promotion; this launcher
# validates the immutable JSON receipt and byte count without rehashing 10 GB
# on the compute node.
#SBATCH --job-name=acs27-source-audit
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/source_audit_extract27_20260919/slurm_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/source_audit_extract27_20260919/slurm_%j.log
set -euo pipefail

ROOT=${KLEVEN_ROOT:-/scratch/td2248/projects/kleven_acs_pilot_20260917}
IN=$ROOT/inputs/ACS/local_extract27_20260919
OUT=$ROOT/output/kleven_acs_pilot/source_audit_extract27_20260919
AUD=$ROOT/forensics/source_audit_extract27_20260919
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
Rscript --vanilla -e 'pkgs <- c("haven","data.table","jsonlite","tidyselect"); stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE))); cat("R_RUNTIME_SMOKE_OK\n")'

step=source_receipt
receipt="$AUD/source_digest_receipt.json"
python3 - "$receipt" "$IN/extract27.dta" <<'PY'
import json, os, re, sys
receipt_path, source_path = sys.argv[1:]
expected = "edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e"
expected_bytes = 9919999546
with open(receipt_path, encoding="utf-8") as handle:
    receipt = json.load(handle)
actual_hash = receipt.get("sha256")
if not isinstance(actual_hash, str) or re.fullmatch(r"[0-9a-f]{64}", actual_hash) is None:
    raise SystemExit("receipt sha256 is not exactly 64 lowercase hex characters")
if actual_hash != expected or receipt.get("expected_sha256") != expected:
    raise SystemExit("receipt sha256 does not match the independently verified source digest")
if int(receipt.get("bytes", -1)) != expected_bytes or os.path.getsize(source_path) != expected_bytes:
    raise SystemExit("source byte count does not match the verified receipt")
print(f"source_sha256={actual_hash} source_size={expected_bytes} source_verification=receipt_and_size")
PY

step=preflight_write
printf 'source_sha256=%s\nsource_size=%s\nsource_verification=receipt_and_size\n' \
  'edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e' \
  '9919999546' > "$OUT/preflight.txt"

step=audit_driver
Rscript --vanilla "$AUD/audit_extract27.R"
step=complete_write
printf 'complete\n' >> "$OUT/preflight.txt"
