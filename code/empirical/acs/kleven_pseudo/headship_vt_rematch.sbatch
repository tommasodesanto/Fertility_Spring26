#!/usr/bin/env bash
# Validation: true restricted rematch of one state (design B) vs filter-and-reweight.
#SBATCH --job-name=acs-headship-vt
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/td2248/projects/kleven_acs_pilot_20260917/output/national_acs_headship_20261112/logs/vt_rematch_%j.log
set -euo pipefail
ROOT=/scratch/td2248/projects/kleven_acs_pilot_20260917
SNAP=$ROOT/snapshots/national_acs_f81814dd
V=$ROOT/output/national_acs_headship_20261112/vt_validation
export HH_GQ="${HH_GQ:?}" STATEFIP=50 KEY_DIR="$ROOT/output/national_acs_headship_20261112/key_local_mac/key"
module load r/4.5.1
export R_LIBS_USER="$ROOT/Rlibs:/scratch/td2248/R_libs" OMP_NUM_THREADS=1
PARTITION_OUTDIR=$V/partitions Rscript --vanilla $ROOT/code_headship_20261112/headship_vt_rematch_prepare.R
( cd $ROOT && PROJECT_ROOT=$SNAP VENDOR_DIR=$SNAP/vendor SOURCE_CONTRACT=$SNAP/code/empirical/acs/kleven_pseudo/source_contract.json \
  HOUSING_RAW=$ROOT/inputs/ACS/local_extract27_20260919/extract27.dta PHASE=smoke STATEFIP=50 \
  PARTITION_OUTDIR=$V/partitions OUTDIR=$V/rematch Rscript --vanilla $SNAP/code/empirical/acs/kleven_pseudo/run_national_acs_match_housing.R )
SAVE_MATCH_TABLE=1 OUTDIR=$V/filter Rscript --vanilla $ROOT/code_headship_20261112/headship_filter_reweight.R
REMATCH_PANEL=$V/rematch/statefip_50/cps_acs_pseudo-panel_housing.rds FILTER_TABLE=$V/filter/states/statefip_50/match_table_B.rds \
  COMPARE_OUT=$V/compare_B.json Rscript --vanilla $ROOT/code_headship_20261112/headship_vt_compare.R
