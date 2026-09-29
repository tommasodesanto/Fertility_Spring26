#!/usr/bin/env bash
# Submit only after lead review and staging: sbatch launch.sh.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:40:00
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
stage=/scratch/td2248/projects/fixed_reference_credit_ge_20260929
old_credit=/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1/credit
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
source="$stage/sources_solve_v1"
results="$stage/ge_results"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
test -d "$reference" && test -d "$source" && test -d "$old_credit"
test ! -e "$results/solve_v1"
mkdir -p "$results"
python - "$source" "$old_credit" "$results/plan_ge.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
source,seed,out=map(Path,sys.argv[1:])
p=json.loads((source/'plan_ge_template.json').read_text())
names=['run_ge.py','natural_credit.py','run_credit.py','run_fixed_price.py',
       'preflight_ge.py','preflight_ge_receipt.json','launch.sh','plan_ge_template.json']
files=[]
for name in names:
    path=source/name
    if not path.is_file(): raise SystemExit('Missing staged GE source/receipt: '+str(path))
    files.append({'path':'/work/ge_source/'+name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
for name in ['receipt.json','conditional_cohort_state.pkl.gz','target_fit.csv','parameters.csv']:
    path=seed/name
    if not path.is_file(): raise SystemExit('Missing credit seed: '+str(path))
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    files.append({'path':'/work/credit_seed/'+name,'sha256':h.hexdigest()})
expected={'receipt.json':p['credit_seed_receipt_sha256'],
          'conditional_cohort_state.pkl.gz':p['credit_seed_checkpoint_sha256']}
for row in files:
    name=Path(row['path']).name
    if row['path'].startswith('/work/credit_seed/') and name in expected and row['sha256']!=expected[name]:
        raise SystemExit('Credit seed SHA differs: '+name)
p['driver_sha256']=files[0]['sha256']
p['files']=files
out.write_text(json.dumps(p,indent=2,sort_keys=True)+'\n')
PY
exec apptainer exec --bind "$reference:$original:ro,$source:/work/ge_source:ro,$results:/work/ge_results:rw,$old_credit:/work/credit_seed:ro" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/ge_source/run_ge.py --plan /work/ge_results/plan_ge.json --output /work/ge_results/solve_v1
