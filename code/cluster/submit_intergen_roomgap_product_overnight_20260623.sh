#!/bin/bash
# Overnight diagnostic split for the June 2026 one-market intergen model.
#
# Submit from the Torch scratch cluster copy:
#   cd $SCRATCH/projects/Fertility_Spring26_20260617_fast/code/cluster
#   bash submit_intergen_roomgap_product_overnight_20260623.sh
#
# Design:
# - 15 concurrent global-DE tasks.
# - 15 concurrent local/random-panel tasks.
# - 120 array tasks per algorithm, so each half can roll through roughly
#   eight 2-hour waves if scheduling is continuous.
# - This is diagnostic only. It does not change model code or target systems.

set -euo pipefail

TARGET_SET="${INTERGEN_TARGET_SET:-candidate_replacement_young_old_roomgap_v1}"
DATE_TAG="${INTERGEN_DATE_TAG:-20260623}"
N_TASKS="${INTERGEN_OVERNIGHT_TASKS:-120}"
CONCURRENCY="${INTERGEN_OVERNIGHT_CONCURRENCY:-15}"
MINUTES="${INTERGEN_MINUTES:-115}"
J="${INTERGEN_J:-16}"
NB="${INTERGEN_NB:-60}"
INCOME_STATES="${INTERGEN_INCOME_STATES:-5}"
N_HOUSE="${INTERGEN_N_HOUSE:-6}"
MAX_ITER_EQ="${INTERGEN_MAX_ITER_EQ:-10}"

GLOBAL_TAG="${INTERGEN_GLOBAL_RUN_TAG:-intergen_roomgap_product_globalde_${DATE_TAG}}"
PANEL_TAG="${INTERGEN_PANEL_RUN_TAG:-intergen_roomgap_product_randompanel_${DATE_TAG}}"

echo "Submitting intergen roomgap/product overnight diagnostic"
echo "target_set=${TARGET_SET}"
echo "arrays=1-${N_TASKS}%${CONCURRENCY} for each algorithm"
echo "J=${J} Nb=${NB} income_states=${INCOME_STATES} n_house=${N_HOUSE} max_iter_eq=${MAX_ITER_EQ}"
echo "minutes_per_task=${MINUTES}"
echo "global_tag=${GLOBAL_TAG}"
echo "panel_tag=${PANEL_TAG}"

GLOBAL_JOB=$(
  INTERGEN_RUN_TAG="${GLOBAL_TAG}" \
  INTERGEN_TARGET_SET="${TARGET_SET}" \
  INTERGEN_GLOBAL_EVALS_PER_TASK="${INTERGEN_GLOBAL_EVALS_PER_TASK:-100000}" \
  INTERGEN_GLOBAL_POP_SIZE="${INTERGEN_GLOBAL_POP_SIZE:-22}" \
  INTERGEN_GLOBAL_MUTATION="${INTERGEN_GLOBAL_MUTATION:-0.85}" \
  INTERGEN_GLOBAL_CROSSOVER="${INTERGEN_GLOBAL_CROSSOVER:-0.70}" \
  INTERGEN_SEED_BASE="${INTERGEN_GLOBAL_SEED_BASE:-2026062300}" \
  INTERGEN_J="${J}" \
  INTERGEN_NB="${NB}" \
  INTERGEN_INCOME_STATES="${INCOME_STATES}" \
  INTERGEN_N_HOUSE="${N_HOUSE}" \
  INTERGEN_MAX_ITER_EQ="${MAX_ITER_EQ}" \
  INTERGEN_MINUTES="${MINUTES}" \
  sbatch --parsable --array="1-${N_TASKS}%${CONCURRENCY}" submit_intergen_housing_fertility_global_de.sh
)

PANEL_JOB=$(
  INTERGEN_RUN_TAG="${PANEL_TAG}" \
  INTERGEN_TARGET_SET="${TARGET_SET}" \
  INTERGEN_CASES_PER_TASK="${INTERGEN_CASES_PER_TASK:-100000}" \
  INTERGEN_WORKERS="${INTERGEN_WORKERS:-1}" \
  INTERGEN_DIAGNOSTIC_BEST="${INTERGEN_DIAGNOSTIC_BEST:-0}" \
  INTERGEN_SEED_BASE="${INTERGEN_PANEL_SEED_BASE:-2026067300}" \
  INTERGEN_J="${J}" \
  INTERGEN_NB="${NB}" \
  INTERGEN_INCOME_STATES="${INCOME_STATES}" \
  INTERGEN_N_HOUSE="${N_HOUSE}" \
  INTERGEN_MAX_ITER_EQ="${MAX_ITER_EQ}" \
  INTERGEN_MINUTES="${MINUTES}" \
  sbatch --parsable --array="1-${N_TASKS}%${CONCURRENCY}" submit_intergen_housing_fertility_twohour_panel.sh
)

echo "Submitted global-DE job: ${GLOBAL_JOB}"
echo "Submitted random-panel job: ${PANEL_JOB}"
echo
echo "Monitor:"
echo "  squeue -j ${GLOBAL_JOB},${PANEL_JOB}"
echo
echo "Expected result roots:"
echo "  results_intergen_housing_fertility_${GLOBAL_TAG}/"
echo "  results_intergen_housing_fertility_${PANEL_TAG}/"
