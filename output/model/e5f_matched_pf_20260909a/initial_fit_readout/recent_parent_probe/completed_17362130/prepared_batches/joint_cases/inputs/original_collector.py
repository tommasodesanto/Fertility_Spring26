"""Read-only 23-case joint-array collector plus explicit smoke reuse. No submits/solves."""
import argparse,base64,csv,datetime,hashlib,importlib.util,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent
DEFAULT_ROUND=pathlib.Path('/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/contracts/joint_round_01')
def read(p):return json.loads(p.read_text())
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def save(p,data):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def csv_write(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def remote(rounddir,pin,jobid):
 # The only launch-gate action invoked by this collector is check-smoke.
 checked=subprocess.run([sys.executable,'-B',str(rounddir/'launch_gate.py'),'check-smoke','--plan-sha256',pin],capture_output=True,text=True)
 if checked.returncode:
  print(json.dumps({'status':'smoke_gate_failed','returncode':checked.returncode,'stdout':checked.stdout,'stderr':checked.stderr,'submitted':False}));return 1
 gate=json.loads(checked.stdout);meta=read(rounddir/'preparation_metadata.json');cases=read(rounddir/'cases.json');out=pathlib.Path(meta['output_root'])
 assert gate['source_pins_verified']==634 and gate['exact_early_equality'] is True
 assert gate['plan_sha256']==pin and len(cases)==24 and len({c['contract_sha256'] for c in cases})==24
 assert gate['proposal_reuse']['case_index']==21
 spec=importlib.util.spec_from_file_location('panel_validator',rounddir/'inputs/panel_validator.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 result={'status':'partial','checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'plan_sha256':pin,'cases_file_sha256':digest(rounddir/'cases.json'),'source_pins_verified':634,'seed_sha256':meta['seed']['sha256'],'smoke_gate_receipt':gate,'array_job_id':jobid,'states':[],'verified_cases':[],'actual_SMM_weights':None,'calibrated_smm':False}
 for case in cases:
  reused=case['index']==21;dest=out/'smoke' if reused else out/'cases'/case['case_id'];rep=dest/('repetition_02' if reused else 'repetition_01')
  state={'index':case['index'],'case_id':case['case_id'],'proposal_contract_sha256':case['contract_sha256'],'smoke_reuse':reused,'remote_output':str(dest)}
  if not (dest/'summary.json').exists():
   state['status']='failed' if any((dest/f).exists() for f in ('failure.json','timeout.json')) else 'pending_or_running'
   for f in ('heartbeat.json','failure.json','timeout.json'):
    if (dest/f).exists():state[f]=read(dest/f)
   result['states'].append(state);continue
  try:
   proposal=read(rounddir/case['contract']);assert digest(rounddir/case['contract'])==case['contract_sha256']
   actual=read(rounddir/'smoke_contract.json') if reused else proposal
   actual_sha=meta['smoke_contract_sha256'] if reused else case['contract_sha256']
   assert read(dest/'contract.json')==dict(actual,contract_sha256=actual_sha,case='new_balanced')
   seed=read(dest/'seed_mapping.json');assert seed['source_files_verified']==634 and seed['initial_psi']==proposal['initial_psi']
   assert actual['source_sha256']==proposal['source_sha256']==meta['source_pins'] and len(proposal['source_sha256'])==634
   assert actual['structural_candidate']==proposal['structural_candidate']
   if reused:
    for key in meta['proposal_reuse']['same_bindings_required']:assert actual[key]==proposal[key]
   top=read(dest/'summary.json');final=read(rep/'summary.json');early=read(rep/'early_measurement.json');ges=read(rep/'stationary_solves.json')
   wrapper=dict(top,repetitions=1,stationary_solves=len(ges),final=final) if reused else top
   params=v.validate_summary(wrapper,proposal,early,list(csv.DictReader((rep/'parameters.csv').open())),ges)
   market=v.validate_market(read(rep/'market_quantity_units.json'))
   checkpoint=rep/'initial_state.pkl.gz';cp_sha=digest(checkpoint);assert cp_sha==final['checkpoint_sha256']
   graphs=sorted((rep/'standard_diagnostics').glob('*.png'));assert [p.name for p in graphs]==sorted(meta['expected_graph_filenames']) and len(graphs)==17
   files={}
   selected=[p for p in sorted(rep.rglob('*')) if p.is_file() and p.suffix in ('.json','.csv')]
   selected += [dest/name for name in ('summary.json','contract.json','seed_mapping.json')]
   for p in selected:
    content=p.read_bytes();files[str(p.relative_to(dest))]={'sha256':hashlib.sha256(content).hexdigest(),'base64':base64.b64encode(content).decode()}
   record={**state,'run_plan_sha256':pin,'input_fingerprint':meta['run_input_fingerprint'],'status':'verified','actual_executed_contract_sha256':actual_sha,'source_commit':'7e872053','proposal':case['proposal'],'parameters':params,'summary':top,'repetition_summary':final,'early_measurement':early,'stationary_solves':len(ges),'elapsed_seconds':top['elapsed_seconds'],'runtime_scope':'entire two-repetition smoke; do not count twice' if reused else 'one fresh case','market_residual':market,'checkpoint':{'path':str(checkpoint),'sha256':cp_sha,'bytes':checkpoint.stat().st_size},'original_graphs':[{'path':str(p),'filename':p.name,'sha256':digest(p),'bytes':p.stat().st_size} for p in graphs],'files':files,'numeric_gates_verified':True}
   result['verified_cases'].append(record);state['status']='verified'
  except Exception as exc:state.update(status='validation_failed',error_type=type(exc).__name__,error=str(exc))
  result['states'].append(state)
 result['status']='complete' if len(result['verified_cases'])==24 else 'partial'
 result['slurm']=subprocess.check_output(['sacct','-j',str(jobid)+',17360699','--format=JobID,State,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,ExitCode,Start,End','-n','-P'],text=True) if jobid else None
 print(json.dumps(result,allow_nan=False));return 0

def ingest(path):
 data=read(path);assert data['status'] in ('partial','complete')
 assert data['plan_sha256']==digest(ROOT.parent/'run_plan.json'),'local/remote run-plan mismatch'
 cases=read(ROOT.parent/'cases.json');assert data['cases_file_sha256']==digest(ROOT.parent/'cases.json')
 assert data['source_pins_verified']==634 and data['smoke_gate_receipt']['exact_early_equality'] is True
 byid={c['case_id']:c for c in cases};collected=ROOT/'array_collected';collected.mkdir(exist_ok=True)
 for item in data['verified_cases']:
  assert item['proposal_contract_sha256']==byid[item['case_id']]['contract_sha256'];dest=collected/item['case_id']
  for name,file in item['files'].items():
   relative=pathlib.Path(name);assert not relative.is_absolute() and '..' not in relative.parts
   content=base64.b64decode(file['base64']);assert hashlib.sha256(content).hexdigest()==file['sha256']
   p=dest/relative;p.parent.mkdir(parents=True,exist_ok=True)
   if p.exists():assert p.read_bytes()==content
   else:p.write_bytes(content)
  receipt={k:v for k,v in item.items() if k!='files'};receipt['file_sha256']={k:v['sha256'] for k,v in item['files'].items()};save(dest/'collection_receipt.json',receipt)
 rows=[]
 for case in cases:
  p=collected/case['case_id']/'collection_receipt.json'
  if p.exists():
   item=read(p);assert item['run_plan_sha256']==data['plan_sha256'],'mixed run-plan fingerprints';rows.append(item)
 raw={'status':'complete' if len(rows)==24 else 'partial','case_count':len(rows),'cases':rows,'original_proposals':24,'fresh_cases':23,'smoke_reuse_index':21,'objective':None,'actual_SMM_weights':None,'calibrated_smm':False,'measurement_note':'Both CPS projections and unavailable moments preserved exactly; no selection or objective.'}
 save(ROOT/'raw_case_summary.json',raw)
 if rows:
  flat=[]
  for row in rows:
   x={k:row[k] for k in ('index','case_id','smoke_reuse','proposal_contract_sha256','actual_executed_contract_sha256','stationary_solves','elapsed_seconds','runtime_scope','market_residual')}
   x.update({'parameter.'+k:v for k,v in row['parameters'].items()})
   x.update({'normalization.'+k:v for k,v in row['repetition_summary']['normalization'].items()})
   for projection,observer in row['early_measurement']['fertility'].items():x.update({'fertility.'+projection+'.'+k:v for k,v in observer['moments'].items()})
   x.update({'housing_wealth.'+k:v for k,v in row['early_measurement']['housing_wealth']['moments'].items()});flat.append(x)
  csv_write(ROOT/'raw_case_summary.csv',flat)
 receipt={k:v for k,v in data.items() if k!='verified_cases'};receipt['case_count_verified']=len(rows);receipt['raw_case_summary_sha256']=digest(ROOT/'raw_case_summary.json');receipt['collector_sha256']=digest(pathlib.Path(__file__));save(ROOT/'array_collection_receipt.json',receipt)
 print(json.dumps({'status':raw['status'],'case_count':len(rows),'raw_case_summary':str(ROOT/'raw_case_summary.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--remote',action='store_true');p.add_argument('--round-dir',type=pathlib.Path,default=DEFAULT_ROUND);p.add_argument('--plan-sha256');p.add_argument('--array-job-id');p.add_argument('--ingest',type=pathlib.Path);a=p.parse_args()
 if a.remote:
  if not a.plan_sha256:raise SystemExit('--plan-sha256 required')
  sys.exit(remote(a.round_dir,a.plan_sha256,a.array_job_id))
 elif a.ingest:ingest(a.ingest)
 else:p.error('choose --remote or --ingest')
