#!/usr/bin/env python3
"""Pinned mechanical run gate. Preparing this file never submits or solves."""
import argparse, csv, datetime, fcntl, hashlib, importlib.util, json, os
from pathlib import Path
import subprocess, sys, time
ROOT=Path(__file__).resolve().parent

def read(p): return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''): h.update(chunk)
 return h.hexdigest()
def save(p,data):
 p.parent.mkdir(parents=True,exist_ok=True)
 temp=p.with_name(p.name+'.'+str(os.getpid())+'.tmp')
 temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n'); os.replace(temp,p)
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def require(ok,message):
 if not ok: raise RuntimeError(message)
def verify_bundle(pin,remote=False):
 require(sha(ROOT/'run_plan.json')==pin,'run-plan fingerprint mismatch')
 plan=read(ROOT/'run_plan.json')
 for name,digest in plan['artifact_sha256'].items(): require(sha(ROOT/name)==digest,'input/artifact fingerprint mismatch: '+name)
 meta=read(ROOT/'preparation_metadata.json'); cases=read(ROOT/'cases.json')
 prototype=read(ROOT/'inputs/panel_smoke_prototype.json')
 require(len(cases)==24 and len({c['contract_sha256'] for c in cases})==24,'24 unique contracts required')
 require([c['index'] for c in cases]==list(range(24)),'case order mismatch')
 analysis=read(ROOT/'inputs/local_analysis.json')
 selected=[(i,p) for i,p in enumerate(analysis['proposals']) if p['ridge_strength'] in (.1,1.) and p['fraction'] in (.5,1.)]
 require(len(selected)==24,'proposal filter mismatch')
 for entry,(index,proposal) in zip(cases,selected):
  require(entry['analysis_index']==index and entry['proposal']==proposal,'analysis order/proposal mismatch')
  require(sha(ROOT/entry['contract'])==entry['contract_sha256'],'case contract fingerprint mismatch')
  c=read(ROOT/entry['contract'])
  require(c['structural_candidate']==proposal['structural_candidate'] and len(c['structural_candidate'])==9,'candidate mismatch')
  require(c['repetitions']==1 and c['maximum_GE_solves']==8,'case loop budget mismatch')
 for name in [c['contract'] for c in cases]+['smoke_contract.json']:
  c=read(ROOT/name)
  require(c['source_sha256']==prototype['source_sha256']==meta['source_pins'] and len(c['source_sha256'])==634,'source pin set mismatch')
  require(c['normalized_checkpoint']==meta['seed']['path'] and c['normalized_checkpoint_sha256']==meta['seed']['sha256'],'seed pin mismatch')
  require(c['source_commit']=='7e872053' and c['normalize'] is True and c['observe_early'] is True,'loop contract mismatch')
  require(c['fertility_normalization']==2.1 and c['initial_psi']==prototype['initial_psi'],'normalization mismatch')
  require(c['maximum_stationary_solves_per_repetition']==8 and c['seconds']==1800,'solve/time budget mismatch')
  require(c['run_input_fingerprint']==meta['run_input_fingerprint'] and 'perturbation' not in c,'input/metadata mismatch')
 smoke=read(ROOT/'smoke_contract.json')
 require(smoke['repetitions']==2 and smoke['maximum_GE_solves']==16,'smoke must use two fresh loops')
 require(smoke['structural_candidate']==cases[21]['proposal']['structural_candidate'] and smoke['proposal_source_index']==47,'smoke selection mismatch')
 require(sha(ROOT/'smoke_contract.json')==meta['smoke_contract_sha256'],'smoke fingerprint mismatch')
 for row in csv.DictReader((ROOT/'inputs/target_fit.csv').open()):
  require(row['actual_weight']=='' and row['loss_contribution']=='','actual SMM weights/contributions must remain unavailable')
 require(read(ROOT/'inputs/raw_case_summary.json')['case_count']==19,'complete panel required')
 if remote:
  snapshot=Path(meta['snapshot'])
  require(ROOT==Path(meta['round_dir']),'wrong remote staging location')
  for name,digest in meta['source_pins'].items(): require(sha(snapshot/name)==digest,'source fingerprint mismatch: '+name)
  require({str(p.relative_to(snapshot)) for p in (snapshot/'code/model').rglob('*.py')}.issubset(meta['source_pins']),'unmanifested model Python source')
  require(sha(Path(meta['seed']['path']))==meta['seed']['sha256'],'seed fingerprint mismatch')
 return meta,cases

def smoke_gate(meta,pin):
 dest=Path(meta['output_root'])/'smoke'; contract=read(ROOT/'smoke_contract.json')
 top=read(dest/'summary.json')
 require(top['status']=='passed_initial_candidate_loop' and top['repetitions']==2 and top['normalized'] is True,'smoke not completed successfully twice')
 require(top['stationary_solves']<=16 and top['elapsed_seconds']<=1800,'smoke exceeded contract budget')
 require(read(dest/'contract.json')==dict(contract,contract_sha256=meta['smoke_contract_sha256'],case='new_balanced'),'actual smoke contract mismatch')
 seed=read(dest/'seed_mapping.json');require(seed['source_files_verified']==634 and seed['initial_psi']==contract['initial_psi'],'smoke seed/source mapping mismatch')
 spec=importlib.util.spec_from_file_location('panel_validator',ROOT/'inputs/panel_validator.py'); validator=importlib.util.module_from_spec(spec);spec.loader.exec_module(validator)
 early=[];reps=[];solves=0;artifacts={}
 for number in (1,2):
  p=dest/f'repetition_{number:02d}'; final=read(p/'summary.json'); observation=read(p/'early_measurement.json');ges=read(p/'stationary_solves.json')
  wrapper=dict(top,repetitions=1,stationary_solves=len(ges),final=final)
  params=validator.validate_summary(wrapper,contract,observation,list(csv.DictReader((p/'parameters.csv').open())),ges)
  validator.validate_market(read(p/'market_quantity_units.json'))
  require(sha(p/'initial_state.pkl.gz')==final['checkpoint_sha256'],'smoke checkpoint fingerprint mismatch')
  early.append(observation);reps.append(final);solves+=len(ges)
  for name in ('summary.json','early_measurement.json','stationary_solves.json','parameters.csv','market_quantity_units.json','initial_state.pkl.gz'):
   artifacts[str((p/name).relative_to(dest))]=sha(p/name)
 require(early[0]==early[1],'fresh smoke early measurements are not exactly equal')
 for name in ('price','normalization','legacy_stationary_moments'):
  require(reps[0][name]==reps[1][name],'fresh smoke repetitions differ: '+name)
 require(top['final']==reps[1] and solves==top['stationary_solves'],'smoke aggregate mismatch')
 graphs=sorted((dest/'repetition_02/standard_diagnostics').glob('*.png'));require(len(graphs)==17,'stable 17-graph packet missing')
 require([p.name for p in graphs]==sorted(meta['expected_graph_filenames']),'standard graph names changed')
 for p in graphs: artifacts[str(p.relative_to(dest))]=sha(p)
 for name in ('summary.json','contract.json','seed_mapping.json'):artifacts[name]=sha(dest/name)
 receipt={'status':'verified','checked_at_utc':now(),'plan_sha256':pin,'input_fingerprint':meta['run_input_fingerprint'],'contract_sha256':meta['smoke_contract_sha256'],'source_pins_verified':634,'seed_sha256':meta['seed']['sha256'],'exact_early_equality':True,'repetitions':2,'stationary_solves':solves,'graph_count':17,'artifact_sha256':artifacts}
 reused=read(ROOT/'contract_21.json')
 for key in ('structural_candidate','source_sha256','normalized_checkpoint','normalized_checkpoint_sha256','initial_psi','normalize','observe_early','fertility_normalization','payroll_tax','housing_supply_elasticity','run_input_fingerprint'):
  require(reused[key]==contract[key],'smoke reuse binding differs: '+key)
 receipt['proposal_reuse']={'case_index':21,'case_id':reused['case_id'],'proposal_contract_sha256':sha(ROOT/'contract_21.json'),'actual_executed_contract_sha256':meta['smoke_contract_sha256'],'output_path':str(dest/'repetition_02'),'original_output_contract_preserved':True,'same_candidate_source_seed_normalization_observers_inputs':True}
 save(Path(meta['output_root'])/'control/smoke_gate_receipt.json',receipt)
 return receipt

def verify_smoke_receipt(meta,pin):
 out=Path(meta['output_root']); receipt=read(out/'control/smoke_gate_receipt.json')
 require(receipt['status']=='verified' and receipt['plan_sha256']==pin and receipt['input_fingerprint']==meta['run_input_fingerprint'] and receipt['exact_early_equality'],'smoke gate receipt mismatch')
 for name,digest in receipt['artifact_sha256'].items():
  # Checkpoints were fully hashed at the submission gate; light files/graphs are rechecked per case.
  if not name.endswith('.pkl.gz'):require(sha(out/'smoke'/name)==digest,'smoke artifact changed: '+name)
 approval=read(out/'control/array_submission.json')
 require(approval['lead_reviewed'] is True and approval['plan_sha256']==pin,'lead review receipt missing')


def submit(kind,meta,pin,lead_reviewed):
 out=Path(meta['output_root']);control=out/'control';control.mkdir(parents=True,exist_ok=True)
 if kind=='array':require(lead_reviewed,'array requires explicit --lead-reviewed');smoke_gate(meta,pin)
 # A durable reservation precedes sbatch. An uncertain response cannot be retried automatically.
 receipt_path=control/(kind+'_submission.json'); reservation=control/(kind+'_submission_reserved.json')
 if receipt_path.exists():
  receipt=read(receipt_path);require(receipt['plan_sha256']==pin,'existing submission fingerprint mismatch');print(json.dumps(receipt));return
 with reservation.open('x') as f:json.dump({'reserved_at_utc':now(),'plan_sha256':pin},f)
 queue=subprocess.check_output(['squeue','-h','-u',os.environ.get('USER','td2248'),'-o','%i|%j'],text=True)
 names={'e5f_joint_smoke','e5f_joint_round01'}
 require(not any(line.split('|')[-1] in names for line in queue.splitlines()),'duplicate round job queued or running')
 receipt={'plan_sha256':pin,'lead_reviewed':bool(lead_reviewed),'requested_at_utc':now(),'original_submission_request':{'cpus':1,'memory_MB':8192,'time':'00:32:00','partition':'automatic','qos':'automatic','account':'torch_pr_570_general','array':'0-20,22-23%24' if kind=='array' else None}}
 # Before array tasks can start they must find this reviewed reservation receipt.
 if kind=='array':save(receipt_path,dict(receipt,status='submission_in_progress'))
 cmd=['sbatch','--parsable','--export=ALL,ROUND_PLAN_SHA256='+pin,str(ROOT/(kind+'.sh'))]
 result=subprocess.run(cmd,text=True,capture_output=True)
 receipt.update(command=cmd,stdout=result.stdout,stderr=result.stderr,returncode=result.returncode,status='submitted' if result.returncode==0 else 'submission_failed_no_retry')
 if result.returncode==0:receipt['job_id']=result.stdout.strip().split(';')[0]
 save(receipt_path,receipt);print(json.dumps(receipt));require(result.returncode==0,'submission failed; reservation retained, do not retry')


def run(kind,meta,cases,pin):
 out=Path(meta['output_root']); control=out/'control';control.mkdir(parents=True,exist_ok=True)
 if kind=='array':
  try:verify_smoke_receipt(meta,pin)
  except Exception as exc:
   save(control/'input_failure.json',{'at_utc':now(),'error':str(exc),'job_id':os.environ.get('SLURM_JOB_ID')});raise
  index=int(os.environ['SLURM_ARRAY_TASK_ID']);require(index!=21,'proposal 21 reuses verified smoke repetition 02; no new solve allowed')
  entry=cases[index];contract_name=entry['contract'];caseid=entry['case_id'];dest=out/'cases'/caseid
  with (control/'stage.lock').open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX)
   require(not (control/'input_failure.json').exists(),'new starts stopped after input/source failure')
   require(len(list(control.glob('candidate_failure_*.json')))<3,'new starts stopped after three candidate failures')
   t=time.time();stagefile=control/'stage.json'
   stage=read(stagefile) if stagefile.exists() else {'first_case_start_unix':t,'first_case_start_utc':now(),'deadline_unix':t+7200,'plan_sha256':pin}
   require(stage['plan_sha256']==pin,'stage fingerprint mismatch')
   require(t+1800<=stage['deadline_unix'],'insufficient time before two-hour stage deadline; new start skipped')
   save(stagefile,stage)
   with (control/('started_'+caseid+'.json')).open('x') as f:json.dump({'started_at_utc':now(),'plan_sha256':pin,'job_id':os.environ.get('SLURM_JOB_ID')},f)
 else:contract_name='smoke_contract.json';caseid='smoke';dest=out/'smoke'
 require(not dest.exists(),'output already exists; no automatic retry or overwrite')
 contract_path=ROOT/contract_name
 cmd=[sys.executable,str(Path(meta['snapshot'])/'code/model/tools/run_e5f_initial_revision_probe.py'),'--contract',str(contract_path),'--contract-sha256',sha(contract_path),'--case','new_balanced','--output',str(dest)]
 started=time.time();process=subprocess.Popen(cmd,cwd=meta['snapshot'])
 while True:
  try:code=process.wait(timeout=max(.01,min(60,1800-(time.time()-started))));break
  except subprocess.TimeoutExpired:
   save(control/('health_'+caseid+'.json'),{'updated_at_utc':now(),'elapsed_seconds':time.time()-started,'job_id':os.environ.get('SLURM_JOB_ID'),'output':str(dest),'driver_heartbeat':read(dest/'heartbeat.json') if (dest/'heartbeat.json').exists() else None})
   if time.time()-started>=1800:
    process.terminate()
    try:process.wait(timeout=5)
    except subprocess.TimeoutExpired:process.kill();process.wait()
    code=124;break
 record={'case_id':caseid,'finished_at_utc':now(),'elapsed_seconds':time.time()-started,'returncode':code,'output':str(dest),'plan_sha256':pin,'job_id':os.environ.get('SLURM_JOB_ID')}
 if code:
  save(control/('candidate_failure_'+caseid+'.json'),record)
 else:
  require((dest/'summary.json').exists(),'driver succeeded without summary')
  record['summary']=read(dest/'summary.json');save(control/('completed_'+caseid+'.json'),record);save(out/'latest_completed_case.json',record)
  if not (out/'best_so_far.json').exists():save(out/'best_so_far.json',{'status':'awaiting_lead_review','reason':'Two diagnostic profiles and projections; no automatic cross-profile selection','actual_SMM_weights':None,'calibrated_smm':False})
 return code

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['check-local','check-remote','check-smoke','submit-smoke','submit-array','run-smoke','run-array']);p.add_argument('--plan-sha256',required=True);p.add_argument('--lead-reviewed',action='store_true');a=p.parse_args()
 meta=read(ROOT/'preparation_metadata.json')
 try:
  meta,cases=verify_bundle(a.plan_sha256,remote=a.action!='check-local')
 except Exception as exc:
  if a.action=='run-array':save(Path(meta['output_root'])/'control/input_failure.json',{'at_utc':now(),'error':str(exc),'job_id':os.environ.get('SLURM_JOB_ID')})
  raise
 if a.action.startswith('check'):
  result=smoke_gate(meta,a.plan_sha256) if a.action=='check-smoke' else {'status':'verified','source_pins':634,'cases':24,'plan_sha256':a.plan_sha256,'remote_verified':a.action=='check-remote'}
  print(json.dumps(result));return 0
 if a.action.startswith('submit'):submit(a.action[7:],meta,a.plan_sha256,a.lead_reviewed);return 0
 return run(a.action[4:],meta,cases,a.plan_sha256)
if __name__=='__main__':sys.exit(main())
