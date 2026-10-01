#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=01:00:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
legacy=/scratch/td2248/projects/fixed_reference_transition_20260928
readiness=/scratch/td2248/projects/transition_readiness_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
sources="$readiness/source/transition_readiness"
test -f "$sources/execution_pins.json"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$readiness/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
apptainer exec --bind "$stage:$original:ro" \
 --bind "$legacy:$original/output/model/fixed_reference_transition_20260928:ro" \
 --bind "$readiness:$original/output/model/transition_readiness_v1:rw" \
 --bind "$sources:$original/code/model/experiments/transition_readiness:ro" \
 --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif bash -s <<'PAYLOAD'
set -euo pipefail
python=/share/apps/anaconda3/2025.06/bin/python
code=code/model/experiments/transition_readiness
packet=output/model/transition_readiness_v1
"$python" -m unittest discover -s "$code" -p test_readiness.py -v > "$packet/preparation/native_loop_smoke.log" 2>&1
"$python" "$code/smoke_exact_loop.py" --output "$packet/preparation/native_control_flow_smoke" > "$packet/preparation/native_control_flow_smoke.log" 2>&1
"$python" "$code/legacy_changed_psi.py" --native-preflight > "$packet/preparation/native_source_preflight.json"
"$python" "$code/legacy_changed_psi.py" --output "$packet/legacy_changed_psi_run"
PAYLOAD
