#!/usr/bin/env bash
# Prepared only; lead review and explicit sbatch submission are required.
#SBATCH --job-name=e5f_fixed_credit_smoke
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'submit with sbatch after lead review' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 NUMBA_DISABLE_JIT=0 NUMBA_CACHE_DIR="/tmp/fixed_credit_numba_cache_${SLURM_JOB_ID}"
mkdir -p "$NUMBA_CACHE_DIR"
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
packet=output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1
hostroot="$stage"
verification="$hostroot/$packet/verification_v1"
mkdir -p "$verification"
receipt="$verification/smoke_${SLURM_JOB_ID}.json"
[[ ! -e "$receipt" ]] || { echo 'refusing to overwrite smoke receipt' >&2; exit 2; }
cd "$hostroot"
sha256sum -c "$hostroot/$packet/source.sha256"
(cd "$hostroot/$packet" && sha256sum -c artifact.sha256)
python - <<'PY'
import hashlib,json,pathlib
root=pathlib.Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project'); packet=root/'output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1'
for name in ('manifest_before.json','manifest_after.json'):
 d=json.loads((packet/'overlay'/name).read_text())
 for f,h in d['sha256'].items():
  p=(root/'code/model/intergen_eqscale_seq_optimized'/f) if name.endswith('before.json') else packet/'overlay'/f
  assert hashlib.sha256(p.read_bytes()).hexdigest()==h, (name,f)
for f in ('prepare_overlay.py','test_contract.py','run_smoke.sh','source.sha256'):
 assert (packet/f).is_file()
PY
python=(apptainer exec --bind "$stage:$original:ro,$NUMBA_CACHE_DIR:$NUMBA_CACHE_DIR,$verification:/work" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python)
"${python[@]}" "$original/$packet/test_contract.py" | tee "$verification/smoke_${SLURM_JOB_ID}.log"
"${python[@]}" - <<'PY'
import hashlib, json, os, pathlib
root=pathlib.Path(os.environ.get('PROJECT_ROOT','/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26'))
p=root/'output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/overlay/manifest_before.json'
data=json.loads(p.read_text()); actual={n:hashlib.sha256((root/'code/model/intergen_eqscale_seq_optimized'/n).read_bytes()).hexdigest() for n in data['sha256']}
assert actual == data['sha256'], 'source pin mismatch'
pathlib.Path('/work/smoke_'+os.environ['SLURM_JOB_ID']+'.json').write_text(json.dumps({'status':'PASS','source_pins':actual,'model_solves':0,'checkpoint_reads':0})+'\n')
PY
