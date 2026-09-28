#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=06:00:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
mode=$1
contract=output/model/fertility_identification_20260928/contract_v1/contract.json
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$stage/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256="$(sha256sum "$stage/project/$contract" | cut -d' ' -f1)"
if [ "$mode" = smoke ]; then
 args=(--stage smoke --output "$original/output/model/fertility_identification_20260928/smoke_v1")
elif [ "$mode" = run ]; then
 approval="$original/output/model/fertility_identification_20260928/approval_v1.json"
 approval_sha=$(sha256sum "$stage/project/output/model/fertility_identification_20260928/approval_v1.json" | cut -d' ' -f1)
 args=(--stage run --output "$original/output/model/fertility_identification_20260928/run_v1" --approval "$approval" --approval-sha256 "$approval_sha")
else
 exit 2
fi
exec apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python code/model/tools/run_e5f_fertility_identification.py --contract "$original/$contract" "${args[@]}"
