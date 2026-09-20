#!/usr/bin/env bash
# Three fixed-price native-entry cohort solves; pure dry-run by default.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
python3 - "$ROOT" <<'PY'
import ast, hashlib, json, os, re, shlex, subprocess, sys, tempfile, time
from pathlib import Path
root=Path(sys.argv[1]); tag=os.environ.get('TAG','income_grid_cohort_v1'); submit=os.environ.get('SUBMIT','0')=='1'; host=os.environ.get('SSH_HOST','torch')
if not re.fullmatch(r'[A-Za-z0-9._-]+',tag): raise SystemExit('invalid TAG')
e='/scratch/td2248/projects/Fertility_Spring26_specification_20260920/'+tag
runtime='/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/finance_dose_refit_v2'
frozen='/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source_v2'
checkpoint='/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/production/selected_verification/evaluation/raw/repetition_02/initial_state.pkl.gz'
base=root/'output/model/native_financing_diagnostic_20260919'; local=base/'specification_followup'/tag
if submit and local.exists(): raise SystemExit('fresh local root required')
stage=local if submit else Path(tempfile.mkdtemp(prefix='income_grid_dryrun_')); stage.mkdir(parents=True,exist_ok=True)
sha=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,v: Path(p).write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
driver=root/'code/model/tools/run_e5f_income_grid_cohort.py'; ast.parse(driver.read_text())
plan_path=base/'overnight/finance_dose_refit_v2/plan.remote.json'; plan=json.loads(plan_path.read_text())
summary=base/'overnight/final_mechanisms/refit_new_income/summary.json'; candidate=base/'earnings_candidate/candidate.json'
helpers={Path(plan[k+'_path']).name:plan[k+'_sha256'] for k in ('constructor','adapter','controller','overnight_controller','factorial_controller')};helpers.update(plan['factorial_helper_sha256'])
assert len(helpers)==9 and plan['source_root']==frozen
constructor=root/'code/model/tools/build_persistent_transitory_income_candidate.py';assert sha(constructor)==helpers[constructor.name]
inputs={
 'checkpoint':{'path':checkpoint,'copy':'checkpoint.pkl.gz','sha256':'b3491eedcee6250cf94833067646d3e6463496cbf5a64bdcabe7b13bc7e89eb2'},
 'summary':{'path':runtime+'/refit_new_income_production/refit_new_income/summary.json','copy':'baseline_summary.json','sha256':sha(summary)},
 'plan':{'path':runtime+'/plan.json','copy':'old_plan.json','sha256':sha(plan_path)},
 'candidate':{'path':e+'/candidate.json','copy':'candidate.json','sha256':sha(candidate)}}
m={'schema':'income_grid_cohort_launch_v3','created':time.time(),'experiment_root':e,'frozen_root':frozen,'runtime_root':runtime,'driver_sha256':sha(driver),'runtime_helpers':helpers,'inputs':inputs,
 'core_hashes':{'solver.py':'2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da','parameters.py':'c0c1c18500fba069152659eaf588c3c895993cdcecb104cee7d6edca6bfae6a5'},
 'design':{'smoke':['5x3'],'production':['9x3','15x3'],'household_solves':3,'case_seconds':900,'smoke_outer_seconds':1050,'production_outer_seconds':2250,'maximum_case_total_seconds':2700,'allocation_minutes':[20,40],'cpu_threads':1,'memory_gb':24,
 'population':'native conditional entry rebuilt; wealth marginal may change','fixed':['annual income parameters','period mapping','three transitory states','all non-income preferences','prices','wealth grid','fiscal objects'],
 'observed_comparable_seconds':84,'time_estimate':'15/27/45 states imply approximately 84/151/252 seconds before extra diagnostics; caps allow JIT and graph overhead','stop':'any source/population/numerical/reproduction/plot failure stops stage; no retries','interpretation':'cohort diagnostic; no target fit, stationary calibration or adoption'}}
write(stage/'launch_manifest.json',m)
remote_stage=r'''import hashlib,json,pathlib,shutil,sys,subprocess
E=pathlib.Path(sys.argv[1]);M=json.loads((E/'launch_manifest.json').read_text());F=pathlib.Path(M['frozen_root'])/'code/model';R=pathlib.Path(M['runtime_root'])/'code/model/tools'
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
for name,pin in M['core_hashes'].items():assert sha(F/'intergen_eqscale_seq_optimized'/name)==pin,name
for name,pin in M['runtime_helpers'].items():assert sha(R/name)==pin,name
for name,pin in M['inputs'].items():
 assert sha(pin['path'])==pin['sha256'],name
 if pathlib.Path(pin['path'])!=E/pin['copy']:shutil.copy2(pin['path'],E/pin['copy'])
assert sha(E/'run_e5f_income_grid_cohort.py')==M['driver_sha256']
D=E/'code/model';D.mkdir(parents=True)
original={str(p.relative_to(F)):sha(p) for p in F.rglob('*.py') if '__pycache__' not in p.parts}
for rel,pin in original.items():
 dest=D/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(F/rel,dest);assert sha(dest)==pin
for name,pin in M['runtime_helpers'].items():shutil.copy2(R/name,D/'tools'/name);assert sha(D/'tools'/name)==pin
shutil.move(str(E/'run_e5f_income_grid_cohort.py'),str(D/'tools/run_e5f_income_grid_cohort.py'))
actual={str(p.relative_to(D)):sha(p) for p in D.rglob('*.py')}
(E/'source_manifest.json').write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
M['source_manifest_sha256']=sha(E/'source_manifest.json');M['frozen_original_source_files']=original
(E/'launch_manifest.json').write_text(json.dumps(M,indent=2,sort_keys=True)+'\n')
for name in ('logs','numba_cache','results'): (E/name).mkdir()
for p in D.rglob('*'):
 if p.is_file():p.chmod(0o444)
for name in ('checkpoint.pkl.gz','baseline_summary.json','old_plan.json','candidate.json','source_manifest.json','launch_manifest.json'): (E/name).chmod(0o444)
print('staged and verified',E)
'''
(stage/'remote_stage.py').write_text(remote_stage)
ast.parse(remote_stage)
for mode,minutes,outer in [('smoke',20,1050),('production',40,2250)]:
 extra='' if mode=='smoke' else ' --smoke-receipt '+shlex.quote(e+'/results/smoke/launch_manifest.json')
 body=f'''#!/usr/bin/env bash
#SBATCH --job-name=grid_{mode}
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:{minutes}:00
#SBATCH --chdir={e}
#SBATCH --output={e}/logs/%x_%j.out
#SBATCH --error={e}/logs/%x_%j.err
#SBATCH --kill-on-invalid-dep=yes
set -euo pipefail
module load anaconda3/2025.06
unset PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0
export NUMBA_CACHE_DIR={shlex.quote(e+'/numba_cache')} MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
cd {shlex.quote(e)}
timeout --signal=TERM --kill-after=30s {outer}s python -B code/model/tools/run_e5f_income_grid_cohort.py --mode {mode} --checkpoint {shlex.quote(e+'/checkpoint.pkl.gz')} --source-root {shlex.quote(e+'/code/model')} --candidate-json {shlex.quote(e+'/candidate.json')} --baseline-summary {shlex.quote(e+'/baseline_summary.json')} --baseline-summary-sha256 {inputs['summary']['sha256']} --output {shlex.quote(e+'/results/'+mode)}{extra}
'''
 (stage/(mode+'.sbatch')).write_text(body);subprocess.run(['bash','-n',str(stage/(mode+'.sbatch'))],check=True)
if not submit: print('dry-run artifacts:',stage);sys.exit()
ssh=['ssh','-o','BatchMode=yes',host]
subprocess.run(ssh+['test ! -e '+shlex.quote(e)+' && mkdir -p '+shlex.quote(e)],check=True)
subprocess.run(['scp','-q',str(stage/'launch_manifest.json'),str(stage/'remote_stage.py'),str(driver),str(candidate),host+':'+e+'/'],check=True)
subprocess.run(ssh+['python3 '+shlex.quote(e+'/remote_stage.py')+' '+shlex.quote(e)],check=True)
subprocess.run(['scp','-q',host+':'+e+'/launch_manifest.json',host+':'+e+'/source_manifest.json',str(stage)+'/'],check=True)
subprocess.run(['scp','-q',str(stage/'smoke.sbatch'),str(stage/'production.sbatch'),host+':'+e+'/'],check=True)
subprocess.run(ssh+['chmod a-w '+shlex.quote(e+'/smoke.sbatch')+' '+shlex.quote(e+'/production.sbatch')],check=True)
r={'remote_root':e,'jobs':{},'status':'submitting','created':time.time()}
for mode in ('smoke','production'):
 dep='' if mode=='smoke' else ' --dependency=afterok:'+r['jobs']['smoke']
 job=subprocess.run(ssh+['sbatch --parsable'+dep+' '+shlex.quote(e+'/'+mode+'.sbatch')],check=True,text=True,capture_output=True).stdout.strip()
 if not job.isdecimal():raise RuntimeError('unexpected submission response '+job)
 r['jobs'][mode]=job;r['status']='smoke_submitted' if mode=='smoke' else 'submitted_afterok';write(stage/'submission.json',r)
print(json.dumps(r))
PY
