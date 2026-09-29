#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
stage=/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929
old_credit=/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
source_v2="$stage/source_v2"
source_v3="$stage/source_v3"
results="$stage/results_v2"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools:/work/elasticity_source"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
test -f "$results/solve_v2/stopped.json"
test -f "$results/plan_v2.json"
test ! -e "$results/solve_v3"
test ! -e "$results/continuation_plan_v3.json"
python - "$source_v3" "$results" <<'PY'
import hashlib,json,sys,subprocess
from pathlib import Path
source,results=map(Path,sys.argv[1:])
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
old=results/'solve_v2'; launch=json.loads((old/'launch.json').read_text())
latest=json.loads((old/'latest_completed.json').read_text())
stop=json.loads((old/'stopped.json').read_text())
if stop.get('reason')!='Conservative observed wall-time forecast exceeds remaining total budget' or latest['lifecycle_solves']!=4:
    raise SystemExit('Original did not stop after four q0 cases solely for forecast')
job=str(launch['slurm_job'])
status=subprocess.run(['sacct','-j',job,'-n','-P','--format=JobID,State,ExitCode'],
    capture_output=True,text=True,check=True,timeout=20)
exact=[line.split('|') for line in status.stdout.splitlines() if line.split('|')[0]==job]
if len(exact)!=1 or exact[0][1:]!=['COMPLETED','0:0']:
    raise SystemExit('Original Slurm job is not completed cleanly')
contract=dict(schema='block0506_elasticity_forecast_continuation_v3',
    original_launch=launch,original_launch_sha256=sha(old/'launch.json'),
    original_deadline_epoch=launch['deadline_epoch'],
    v2_plan_sha256=sha(results/'plan_v2.json'),
    v2_driver_sha256=sha(source.parent/'source_v2/run_elasticity.py'),
    continuation_driver_sha256=sha(source/'continue_remaining.py'),
    original_job_terminal=dict(job_id=job,state='COMPLETED',exit_code='0:0'),
    completed_records=latest['completed'],no_new_deadline=True,
    maximum_total_lifecycle_solves=12,case_seconds=600,thread_count=1,memory_gib=16)
path=results/'continuation_plan_v3.json'
path.write_text(json.dumps(contract,indent=2,sort_keys=True)+'\n')
path.chmod(0o444)
PY
exec apptainer exec --bind "$reference:$original:ro,$source_v2:/work/elasticity_source:ro,$source_v3:/work/elasticity_continue_source:ro,$results:/work/elasticity_results:rw,$results/solve_v2:/work/elasticity_completed:ro,$old_credit:/work/credit_seed:ro" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/elasticity_continue_source/continue_remaining.py --contract /work/elasticity_results/continuation_plan_v3.json --plan /work/elasticity_results/plan_v2.json --original /work/elasticity_completed --output /work/elasticity_results/solve_v3
