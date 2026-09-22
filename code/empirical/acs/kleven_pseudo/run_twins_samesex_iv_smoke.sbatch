#!/usr/bin/env bash
# Twin1/SameSex2 IV smoke: single small state (VT, statefip 50, 121,537 rows)
# through the exact production loop path before the national submission.
#SBATCH --job-name=acs-twins-samesex-smoke
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_smoke_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_smoke_%j.log
set -euo pipefail
ROOT_REPO=/scratch/td2248/projects/kleven_acs_pilot_20260917
OUTDIR="${OUTDIR:?set OUTDIR to a fresh directory}"
module load r/4.5.1
export R_LIBS_USER="$ROOT_REPO/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT_REPO"
ROOT="$ROOT_REPO/code/empirical/acs/kleven_pseudo" \
PARTITION_DIR="$ROOT_REPO/output/national_acs_source_stage_20260921/partitions" \
STATEFIP_LIST=50 \
OUTDIR="$OUTDIR" \
STATUS_PATH="$OUTDIR/status_heartbeat.json" \
SAMPLE_LABEL=young \
Rscript --vanilla code/empirical/acs/kleven_pseudo/run_twins_samesex_iv.R
