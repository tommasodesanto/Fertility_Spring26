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
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$stage/numba_cache"
export PYTHONPATH="$original/code/model/tools:$original/tmp/e5f_overnight_local_20260927/portable/tools_v4${PYTHONPATH:+:$PYTHONPATH}"
base=output/model/fertility_identification_20260928
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
if [ "$mode" = check ]; then
 exec apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python -c 'from pathlib import Path; import time; import run_e5f_fertility_identification_resume as r; b=Path("output/model/fertility_identification_20260928"); c,o=r.original.verify(b/"contract_v1/contract.json"); x,y=r.authenticate(c,o,r.core.read(b/"resume_manifest_v1.json"),r.core.sha(b/"contract_v1/contract.json"),time.time()); print("Authenticated original completed records:",len(x),"; zero solves")'
fi
[ "$mode" = run ] || exit 2
manifest_sha=$(sha256sum "$stage/project/$base/resume_manifest_v1.json" | cut -d' ' -f1)
approval_sha=$(sha256sum "$stage/project/$base/approval_v1.json" | cut -d' ' -f1)
wrapper_sha=$(sha256sum "$stage/project/code/model/tools/run_e5f_fertility_identification_resume.py" | cut -d' ' -f1)
exec apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python code/model/tools/run_e5f_fertility_identification_resume.py --contract "$original/$base/contract_v1/contract.json" --output "$original/$base/resume_v1" --approval "$original/$base/approval_v1.json" --approval-sha256 "$approval_sha" --resume-manifest "$original/$base/resume_manifest_v1.json" --resume-manifest-sha256 "$manifest_sha" --wrapper-sha256 "$wrapper_sha"
