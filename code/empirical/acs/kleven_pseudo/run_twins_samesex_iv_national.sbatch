#!/usr/bin/env bash
# Twin1/SameSex2 IV national production run: all 50 states + DC through the
# same loop path validated by the smoke job. Checkpoints after each outcome.
#SBATCH --job-name=acs-twins-samesex-national
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_national_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_twins_samesex_national_%j.log
set -euo pipefail
ROOT_REPO=/scratch/td2248/projects/kleven_acs_pilot_20260917
OUTDIR="${OUTDIR:?set OUTDIR to a fresh directory}"
STATEFIP_LIST="1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56"
module load r/4.5.1
export R_LIBS_USER="$ROOT_REPO/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT_REPO"
ROOT="$ROOT_REPO/code/empirical/acs/kleven_pseudo" \
PARTITION_DIR="$ROOT_REPO/output/national_acs_source_stage_20260921/partitions" \
STATEFIP_LIST="$STATEFIP_LIST" \
OUTDIR="$OUTDIR" \
STATUS_PATH="$OUTDIR/status_heartbeat.json" \
SAMPLE_LABEL=young \
Rscript --vanilla code/empirical/acs/kleven_pseudo/run_twins_samesex_iv.R
