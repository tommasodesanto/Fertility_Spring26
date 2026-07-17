#!/bin/bash
# Strictly collect the five conditional theta1 profile cells.
#SBATCH --job-name=ihfstdbqpc
#SBATCH --output=logs/slurm_ihf_stdbqpc_%j.out
#SBATCH --error=logs/slurm_ihf_stdbqpc_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"

RUN_ROOT="$PROJECT_ROOT/output/model/intergen_standard_bequest_recalibration_20260716"
WINNER_JSON="${STANDARD_BEQUEST_WINNER_JSON:-$RUN_ROOT/report/results.json}"
PROFILE_ROOT="${STANDARD_BEQUEST_PROFILE_ROOT:-$RUN_ROOT/theta1_profile}"
OUTDIR="${STANDARD_BEQUEST_PROFILE_REPORT_DIR:-$RUN_ROOT/theta1_profile_report}"
exec "$PYTHON_BIN" "$MODEL_DIR/tools/collect_intergen_standard_bequest_theta1_profile.py" \
  --winner-json "$WINNER_JSON" \
  --results-dir "$PROFILE_ROOT" \
  --outdir "$OUTDIR"
