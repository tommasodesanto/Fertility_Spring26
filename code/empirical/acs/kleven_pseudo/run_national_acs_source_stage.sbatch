#!/usr/bin/env bash
# Inventory/partition smoke only. Full national matching requires lead review of
# source_schema_receipt.json and partition_receipt.json.
#SBATCH --job-name=acs-source-stage
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_source_stage_%j.log
#SBATCH --error=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_source_stage_%j.log
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
OUTDIR="${OUTDIR:?set OUTDIR to a fresh directory}"
STATEFIP_LIST="1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56"
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT"
PHASE=partition_smoke STATEFIP="$STATEFIP_LIST" Rscript --vanilla code/empirical/acs/kleven_pseudo/run_national_acs_source_stage.R
