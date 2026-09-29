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
source_v4="$stage/source_v4"
results="$stage/results_v2"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools:/work/elasticity_source"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
test -f "$results/solve_v2/stopped.json"
test -f "$results/plan_v2.json"
test ! -e "$results/solve_v4"
test ! -e "$results/continuation_plan_v4.json"
python - "$source_v4" "$results" <<'PY'
import hashlib,json,sys,subprocess
from pathlib import Path
source,results=map(Path,sys.argv[1:])
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
old=results/'solve_v2'; launch=json.loads((old/'launch.json').read_text())
if sha(source/'continue_remaining.py')!='a2f66adaa6331ec4217f7e5fe07e64978b6c85fcf3b046642f2d87f5d595d0a5' or sha(source.parent/'source_v2/run_elasticity.py')!='5d461d4999c428c77438ee43dd63db0e69ee0dd10356b0e09d65c62880adc5ca':
    raise SystemExit('Frozen v4 controller or v2 driver source differs')
if str(launch['slurm_job'])!='18801318' or launch['deadline_epoch']!=1790700222.6872504 or sha(results/'plan_v2.json')!='6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f':
    raise SystemExit('Original job, deadline, or plan identity differs')
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
contract=dict(schema='block0506_elasticity_forecast_continuation_v4',
    original_launch=launch,original_launch_sha256=sha(old/'launch.json'),
    original_deadline_epoch=launch['deadline_epoch'],
    v2_plan_sha256=sha(results/'plan_v2.json'),
    v2_driver_sha256=sha(source.parent/'source_v2/run_elasticity.py'),
    continuation_driver_sha256=sha(source/'continue_remaining.py'),
    original_job_terminal=dict(job_id=job,state='COMPLETED',exit_code='0:0'),
    completed_records=latest['completed'],no_new_deadline=True,
    maximum_total_lifecycle_solves=12,case_seconds=600,thread_count=1,memory_gib=16)
path=results/'continuation_plan_v4.json'
path.write_text(json.dumps(contract,indent=2,sort_keys=True)+'\n')
path.chmod(0o444)
PY
exec apptainer exec --bind "$reference:$original:ro,$source_v2:/work/elasticity_source:ro,$source_v4:/work/elasticity_continue_source:ro,$results:/work/elasticity_results:rw,$results/solve_v2:/work/elasticity_completed:ro,$old_credit:/work/credit_seed:ro" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/elasticity_continue_source/continue_remaining.py --contract /work/elasticity_results/continuation_plan_v4.json --plan /work/elasticity_results/plan_v2.json --original /work/elasticity_completed --output /work/elasticity_results/solve_v4
