#!/usr/bin/env bash
# Initial user-authorized launch payload; monitoring actions must remain read-only.
# This script prepares a pinned plan, then delegates execution to the existing launcher.
# It never submits, retries, or cancels a Slurm job.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=13:00:00
# The allocation supplies memory headroom; numerical threads and estimator process stay single-threaded.
set -euo pipefail

ROOT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
overlay=/scratch/td2248/projects/fixed_reference_transition_20260928
relative=output/model/fixed_reference_transition_20260928
packet="$overlay/four_shock_v1/launch_v2"
controller="$packet/source_launch/code/cluster/run_e5f_preference_estimation_batch.py"
launcher="$packet/source_launch/code/cluster/submit_e5f_preference_estimation.sh"
kind="${E5F_PE_KIND:?set E5F_PE_KIND to four_successive or one_permanent}"
controller_sha="${E5F_PE_CONTROLLER_SHA256:?set controller SHA-256}"
launcher_sha="${E5F_PE_LAUNCHER_SHA256:?set launcher SHA-256}"
gated_sha="${E5F_PE_GATED_SHA256:?set gated script SHA-256}"
case "$kind" in four_successive|one_permanent) ;; *) echo 'invalid E5F_PE_KIND' >&2; exit 2 ;; esac

# All three implementation pins are checked before plan creation or environment setup.
test "$(sha256sum "$controller" | awk '{print $1}')" = "$controller_sha"
test "$(sha256sum "$launcher" | awk '{print $1}')" = "$launcher_sha"
test "$(sha256sum "$0" | awk '{print $1}')" = "$gated_sha"

plan_rel="$relative/four_shock_v1/launch_v2/plans/$kind.json"
output_rel="$relative/four_shock_v1/launch_v2/runs/$kind"
plan="$ROOT/$plan_rel"
plan_physical="$packet/plans/$kind.json"
output_physical="$packet/runs/$kind"
source_rel="$relative/four_shock_v1/launch_v2/source_launch/code/model/tools"
logical_packet="$ROOT/$relative/four_shock_v1/launch_v2"
logical_controller="$logical_packet/source_launch/code/cluster/run_e5f_preference_estimation_batch.py"
test ! -e "$plan_physical"
test ! -e "$output_physical"
mkdir -p "$packet/plans"

module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export NUMBA_CACHE_DIR="$overlay/four_shock_v1/numba_cache"
apptainer exec --bind "$stage/project:$ROOT:ro" --bind "$overlay:$ROOT/$relative:rw" --pwd "$ROOT" \
  /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$logical_controller" create-plan \
  --kind "$kind" --source-dir "$ROOT/$source_rel" --readiness "$ROOT/$relative/four_shock_v1/launch_v2/smoke/readiness.json" \
  --empirical-blocks "$ROOT/$relative/four_shock_v1/estimation_inputs/empirical_blocks.csv" \
  --annual "$ROOT/$relative/four_shock_v1/estimation_inputs/annual_fertility_2007_2023.csv" --output "$plan" \
  --reference-manifest-sha256 147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4 \
  --housing fixed_stock --horizon 104 --horizon 128 --seed-horizon 12 --perturbed-date 5 --log-step 1e-5 \
  --cache-max-bytes 68719476736 --jacobian-receipt "$ROOT/$relative/four_shock_v1/launch_v2/smoke/native/jacobian_seed/measured/receipt.json" \
  --fit-evaluations 12 --endpoint-evaluations 24 --path-evaluations 16 \
  --total-seconds 43200 --candidate-seconds 36000 --endpoint-seconds 3600 \
  --mapping-seconds 5400 --path-seconds 21600 --jacobian-seconds 3000 --enable-execution

plan_sha="$(sha256sum "$plan_physical" | awk '{print $1}')"
export E5F_PE_CONTROLLER_RELATIVE="$relative/four_shock_v1/launch_v2/source_launch/code/cluster/run_e5f_preference_estimation_batch.py"
export E5F_PE_SOURCE_RELATIVE="$source_rel"
export E5F_PE_PLAN_RELATIVE="$plan_rel" E5F_PE_OUTPUT_RELATIVE="$output_rel" E5F_PE_PLAN_SHA256="$plan_sha"
export E5F_PE_CONTROLLER_SHA256="$controller_sha" E5F_PE_LAUNCHER_SHA256="$launcher_sha"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
exec bash "$launcher"
