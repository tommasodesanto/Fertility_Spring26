#!/usr/bin/env bash
# Submit only after lead review of replay_supply.py and plan.json.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
ge=/scratch/td2248/projects/fixed_reference_credit_ge_20260929
stage=/scratch/td2248/projects/fixed_reference_supply_20260929
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
test -f "$reference/output/model/fertility_identification_20260928/resume_v1/selected_export/primary/initial_state.pkl.gz"
test -f "$ge/ge_results/solve_v1/selected_repeat/conditional_cohort_state.pkl.gz"
test -f "$stage/sources_v1/replay_supply.py"
test -f "$stage/sources_v1/plan.json"
test ! -e "$stage/results_v1/replay"
mkdir -p "$stage/results_v1"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools:/work/ge_source"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
exec apptainer exec --bind "$reference:$original:ro,$ge/sources_solve_v1:/work/ge_source:ro,$ge/ge_results:/work/ge_results:ro,$stage/sources_v1:/work/supply_source:ro,$stage/results_v1:/work/supply_results:rw" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/supply_source/replay_supply.py --plan /work/supply_source/plan.json --output /work/supply_results/replay
