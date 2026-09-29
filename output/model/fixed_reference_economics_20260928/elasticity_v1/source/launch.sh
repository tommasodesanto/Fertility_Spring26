#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:40:00
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
stage=/scratch/td2248/projects/fixed_reference_elasticity_20260929
old_credit=/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
source="$stage/source_v1"
results="$stage/results_v1"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools:/work/elasticity_source"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
test -d "$reference" && test -d "$source" && test -d "$old_credit"
test -f "$results/preflight_v1/preflight.json"
test ! -e "$results/solve_v1"
test ! -e "$results/plan_v1.json"
python - "$source" "$results" "$old_credit" <<'PY'
import hashlib,json,sys
from pathlib import Path
source,results,old=map(Path,sys.argv[1:])
plan=json.loads((source/'plan_template.json').read_text())
files=[]
for item in source.iterdir():
    if item.is_file():
        files.append({'path':'/work/elasticity_source/'+item.name,
                      'sha256':hashlib.sha256(item.read_bytes()).hexdigest()})
pre=results/'preflight_v1/preflight.json'
receipt=json.loads(pre.read_text())
if receipt['status']!='ready' or receipt['model_solves']!=0: raise SystemExit('Preflight did not pass')
files.append({'path':'/work/elasticity_results/preflight_v1/preflight.json',
              'sha256':hashlib.sha256(pre.read_bytes()).hexdigest()})
for regime in ('grid_control','credit'):
    for name in ('receipt.json','conditional_cohort_state.pkl.gz','target_fit.csv','parameters.csv'):
        f=old/regime/name
        if not f.is_file(): raise SystemExit('Missing old replay input: '+str(f))
        h=hashlib.sha256()
        with f.open('rb') as stream:
            for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
        files.append({'path':'/work/credit_seed/'+regime+'/'+name,'sha256':h.hexdigest()})
plan['driver_sha256']=next(x['sha256'] for x in files if x['path'].endswith('/run_elasticity.py'))
plan['files']=files
(results/'plan_v1.json').write_text(json.dumps(plan,indent=2,sort_keys=True)+'\n')
PY
exec apptainer exec --bind "$reference:$original:ro,$source:/work/elasticity_source:ro,$results:/work/elasticity_results:rw,$old_credit:/work/credit_seed:ro" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/elasticity_source/run_elasticity.py --plan /work/elasticity_results/plan_v1.json --output /work/elasticity_results/solve_v1
