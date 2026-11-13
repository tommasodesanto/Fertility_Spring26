#!/usr/bin/env bash
# Headship variant of the national Twin1/SameSex2 IV run (2026-11-12).
# Same resources and env as run_twins_samesex_iv_national_128g.sbatch.
#SBATCH --job-name=acs-hs-iv
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_headship_20261112/logs/iv_%x_%j.log
set -euo pipefail
ROOT_REPO=/scratch/td2248/projects/kleven_acs_pilot_20260917
CODE_ROOT=$ROOT_REPO/code_headship_20261112
MODE="${HEADSHIP_MODE:?restrict or outcome}"
OUTDIR=$ROOT_REPO/output/national_acs_headship_20261112/iv_$MODE
if [ -d "$OUTDIR" ]; then echo "OUTDIR exists: $OUTDIR" >&2; exit 1; fi
module load r/4.5.1
export R_LIBS_USER="$ROOT_REPO/Rlibs:/scratch/td2248/R_libs" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$ROOT_REPO"
ROOT="$CODE_ROOT" PARTITION_DIR="$ROOT_REPO/output/national_acs_source_stage_20260921/partitions" \
STATEFIP_LIST="1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56" \
OUTDIR="$OUTDIR" STATUS_PATH="$OUTDIR/status_heartbeat.json" SAMPLE_LABEL=young \
HEADSHIP_MODE="$MODE" KEY_DIR="$ROOT_REPO/output/national_acs_headship_20261112/key_local_mac/key" HH_GQ_FILE="$CODE_ROOT/hh_gq.txt" \
Rscript --vanilla "$CODE_ROOT/headship_twins_samesex_iv.R"
