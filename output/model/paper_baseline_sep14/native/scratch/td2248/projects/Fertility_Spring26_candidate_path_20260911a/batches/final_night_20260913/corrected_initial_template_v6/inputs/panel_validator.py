"""Read-only remote panel verifier and local light-artifact ingester. No solves."""
import argparse,base64,csv,datetime,hashlib,io,json,math,pathlib,subprocess
REMOTE=pathlib.Path('/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053')
NAMES=('beta_annual','kappa_fert','kappa_fert_continuation','chi','H0','theta0','theta1','first_birth_fixed_cost','h_P')
BOUNDS={'beta_annual':(.94,.9995),'kappa_fert':(.02,50.),'kappa_fert_continuation':(.02,50.),'chi':(.10,5.),'H0':(.20,80.),'theta0':(0.,8.),'theta1':(.02,16.),'first_birth_fixed_cost':(0.,8.),'h_P':(.10,2.30)}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def finite(x):return isinstance(x,(float,int)) and not isinstance(x,bool) and math.isfinite(x)
def same(a,b):return math.isclose(float(a),float(b),rel_tol=0,abs_tol=1e-12)
def validate_market(quantity):
 q=quantity['aggregate_quantities'];d=q['housing_demand'];s=q['housing_supply']
 assert len(d)==len(s)==1 and all(finite(x) for x in d+s) and s[0]>0
 residual=max(abs(x-y)/max(y,1e-12) for x,y in zip(d,s))
 assert residual<=2e-4
 assert same(q['aggregate_housing_demand'],sum(d)) and same(q['aggregate_housing_supply'],sum(s))
 assert same(q['aggregate_housing_excess'],sum(d)-sum(s))
 assert same(q['aggregate_owner_demand']+q['aggregate_rental_demand'],sum(d))
 return residual
def validate_summary(summary,contract,early,parameter_rows,ges):
 assert summary['status']=='passed_initial_candidate_loop'
 assert summary['normalized'] is True and summary['repetitions']==1
 assert summary['calibrated_smm'] is False and summary['perfect_foresight_solved'] is False
 assert 1<=summary['stationary_solves']<=contract['maximum_stationary_solves_per_repetition']==8
 assert len(ges)==summary['stationary_solves']
 final=summary['final']; norm=final['normalization']; gates=final['operator_gates']
 assert norm['target']==2.1 and finite(norm['completed_fertility']) and abs(norm['completed_fertility']-2.1)<=5e-4
 assert same(norm['absolute_gap'],abs(norm['completed_fertility']-2.1))
 for k in ('stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1','mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):
  assert finite(gates[k]) and abs(gates[k])<=5e-9,k
 assert finite(gates['zero_entry_mass_accounting_residual']) and abs(gates['zero_entry_mass_accounting_residual'])<=2e-8
 assert finite(gates['stationary_feasibility_projection_mass']) and 0<=gates['stationary_feasibility_projection_mass']<=1e-6
 budget=final['household_budget'];assert budget['budget_excess_mass']==0 and budget['maximum_occupied_excess']<=budget['budget_tolerance']<=1e-9
 arrays=final['policy_array_gates'];assert arrays['occupied_negative_steps']==0
 for stats in arrays['probabilities'].values():
  assert stats['nonfinite']==0 and finite(stats['minimum']) and finite(stats['maximum']) and 0<=stats['minimum']<=stats['maximum']<=1
 def fiscal_ok(f):
  assert f['fiscal_gate'] and f['marginal_gate']
  assert finite(f['normalized_age_income_max_gap']) and f['normalized_age_income_max_gap']<=f['marginal_tolerance']<=1e-9
  a=f['actual_accounts'];assert finite(a['scaled_pension_budget_residual']) and abs(a['scaled_pension_budget_residual'])<=f['fiscal_tolerance']<=1e-6
  assert a['payroll_tax_rate']==.179
  assert same(a['pension_budget_residual'],a['payroll_tax_revenue']-a['pension_outlays'])
 fiscal_ok(final['fiscal'])
 for row in ges:
  assert row['status']=='completed' and finite(row['market_error']) and 0<=row['market_error']<=2.5e-5
  fiscal_ok(row['fiscal'])
 params={r['parameter']:float(r['estimate']) for r in parameter_rows}
 assert len(params)==len(parameter_rows)
 assert set(contract['structural_candidate'])==set(NAMES)
 for row in parameter_rows:
  name=row['parameter'];estimate=float(row['estimate']);assert finite(estimate)
  if name in NAMES:
   low,high=BOUNDS[name]
   assert (float(row['lower']),float(row['upper']))==(low,high)
   assert low<=estimate<=high and same(estimate,contract['structural_candidate'][name]),name
   assert (str(row['near_bound']).lower()=='true')==(min(estimate-low,high-estimate)<=.01*(high-low))
 for name,value in {'hbar_child_rooms':0.,'payroll_tax':.179,'housing_supply_elasticity':.63,'alpha_cons':.733,'sigma':2.}.items():assert params[name]==value
 assert same(params['psi_child'],norm['psi_child'])
 assert early==final['early_measurement']
 assert early['calibrated_smm'] is False and early['weights_assigned'] is False and early['empirical_target_contract_activated'] is False
 assert set(early['fertility'])=={'uniform_birth_time','constant_post_cell'}
 for projection,observer in early['fertility'].items():
  assert observer['metadata']['age_projection']==projection and observer['metadata']['production_smm_eligible'] is False
  assert observer['metadata']['weights_or_standard_errors_adopted'] is False
  assert all(finite(v) for v in observer['moments'].values())
  assert abs(sum(observer['parity_shares_40_44'].values())-1)<=2e-12
  assert observer['accounting']['maximum_age_mass_error']<=observer['metadata']['mass_atol']==2e-10
  assert observer['accounting']['maximum_parity_flow_error']<=observer['metadata']['flow_atol']==2e-10
 housing=early['housing_wealth'];assert housing['production_eligible'] is False and housing['target_contract_activated'] is False
 assert housing['age_projection']=='uniform_within_age_cell'
 for row in housing['rows']:
  value=row['model_value'];assert (value is None and row['available'] is False) or (finite(value) and row['available'] is True)
  assert value==housing['moments'][row['moment']]
 return params

def remote(skip):
 cases=read(REMOTE/'contracts/cases.json');assert len(cases)==19
 assert len({c['case_id'] for c in cases})==len({c['contract_sha256'] for c in cases})==19
 assert [c['index'] for c in cases]==list(range(19))
 contracts={}
 for case in cases:
  path=REMOTE/'contracts'/case['contract'];assert digest(path)==case['contract_sha256']
  c=read(path);assert c['case_id']==case['case_id'] and c['source_commit']=='7e872053'
  assert c['normalize'] is True and c['observe_early'] is True and c['repetitions']==1 and c['seconds']==1800
  assert c['perturbation']==case['perturbation'];contracts[case['case_id']]=c
 baseline=contracts['baseline'];source=baseline['source_sha256'];assert len(source)==634
 for c in contracts.values():
  assert c['source_sha256']==source and c['normalized_checkpoint']==baseline['normalized_checkpoint'] and c['normalized_checkpoint_sha256']==baseline['normalized_checkpoint_sha256']
  point=c['structural_candidate'];assert set(point)==set(NAMES)
  if c['case_id']!='baseline':
   perturbation=c['perturbation'];parameter=perturbation['parameter']
   assert parameter in NAMES and same(perturbation['base'],baseline['structural_candidate'][parameter])
   assert same(point[parameter],perturbation['value'])
   assert same(point[parameter]-perturbation['base'],perturbation['actual_change'])
   assert [name for name in NAMES if point[name]!=baseline['structural_candidate'][name]]==[parameter]
   assert (perturbation['actual_change']>0)==c['case_id'].endswith('_plus')
 for name,pin in source.items():assert digest(REMOTE/name)==pin,name
 assert digest(pathlib.Path(baseline['normalized_checkpoint']))==baseline['normalized_checkpoint_sha256']
 assert {str(p.relative_to(REMOTE)) for p in (REMOTE/'code/model').rglob('*.py')}.issubset(source)
 result={'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_commit':'7e872053','source_pins_verified':634,'seed_sha256_verified':baseline['normalized_checkpoint_sha256'],'all_19_unique_contracts_verified':True,'cases_file_sha256':digest(REMOTE/'contracts/cases.json'),'cases':[],'newly_verified':[]}
 for case in cases:
  caseid=case['case_id'];path=REMOTE/'output/panel'/caseid;state={'case_id':caseid,'index':case['index'],'contract_sha256':case['contract_sha256']}
  for name in ('heartbeat.json','failure.json','timeout.json'):
   if (path/name).exists():state[name]=read(path/name)
  if caseid in skip:
   state['status']='previously_verified';result['cases'].append(state);continue
  if not (path/'summary.json').exists():
   state['status']='failed' if 'failure.json' in state or 'timeout.json' in state else 'pending_or_running';result['cases'].append(state);continue
  try:
   contract=contracts[caseid];actual=read(path/'contract.json')
   assert actual=={**contract,'contract_sha256':case['contract_sha256'],'case':'new_balanced'}
   seed=read(path/'seed_mapping.json');assert seed['source_files_verified']==634 and seed['initial_psi']==contract['initial_psi']
   summary=read(path/'summary.json');early=read(path/'repetition_01/early_measurement.json')
   parameters=list(csv.DictReader((path/'repetition_01/parameters.csv').open()))
   ges=read(path/'repetition_01/stationary_solves.json');values=validate_summary(summary,contract,early,parameters,ges)
   validate_market(read(path/'repetition_01/market_quantity_units.json'))
   checkpoint=path/'repetition_01/initial_state.pkl.gz';checkpoint_sha=digest(checkpoint)
   assert checkpoint_sha==summary['final']['checkpoint_sha256']
   graphs=sorted((path/'repetition_01/standard_diagnostics').glob('*.png'));assert len(graphs)==17
   files={}
   for p in sorted(path.rglob('*')):
    if p.is_file() and p.suffix in ('.json','.csv'):
     data=p.read_bytes();files[str(p.relative_to(path))]={'sha256':hashlib.sha256(data).hexdigest(),'bytes_base64':base64.b64encode(data).decode(),'size_bytes':len(data)}
   graph_receipt=[{'path':str(p),'filename':p.name,'sha256':digest(p),'bytes':p.stat().st_size} for p in graphs]
   record={**state,'status':'verified','parameter_values':values,'summary':summary,'early_measurement':early,'checkpoint':{'path':str(checkpoint),'sha256':checkpoint_sha,'bytes':checkpoint.stat().st_size},'original_graphs':graph_receipt,'files':files}
   result['newly_verified'].append(record);state['status']='verified'
  except Exception as exc:
   state.update(status='collection_validation_failed',error_type=type(exc).__name__,error=str(exc))
  result['cases'].append(state)
 result['slurm']=subprocess.check_output(['sacct','-j','17358647','--format=JobID,State,Elapsed,TotalCPU,MaxRSS,ExitCode','-n','-P'],text=True)
 print(json.dumps(result,allow_nan=True))

def ingest(path):
 root=pathlib.Path(__file__).resolve().parent;collected=root/'collected';collected.mkdir(exist_ok=True)
 data=read(path);cases=read(root/'cases.json');byid={c['case_id']:c for c in cases}
 assert data['cases_file_sha256']==digest(root/'cases.json')
 assert data['all_19_unique_contracts_verified'] and data['source_pins_verified']==634
 for case in cases:assert digest(root/case['contract'])==case['contract_sha256']
 for case in data['newly_verified']:
  caseid=case['case_id'];assert case['contract_sha256']==byid[caseid]['contract_sha256']
  dest=collected/caseid;dest.mkdir(exist_ok=True)
  for name,entry in case['files'].items():
   rel=pathlib.Path(name);assert not rel.is_absolute() and '..' not in rel.parts
   content=base64.b64decode(entry['bytes_base64']);assert hashlib.sha256(content).hexdigest()==entry['sha256']
   file=dest/rel;file.parent.mkdir(parents=True,exist_ok=True)
   if file.exists():assert file.read_bytes()==content
   else:file.write_bytes(content)
  receipt={k:v for k,v in case.items() if k!='files'}
  receipt['file_sha256']={k:v['sha256'] for k,v in case['files'].items()}
  (dest/'collection_receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=True)+'\n')
 (collected/'latest_slurm.txt').write_text(data['slurm'])
 rows=[]
 for case in cases:
  file=collected/case['case_id']/'collection_receipt.json'
  if not file.exists():continue
  receipt=read(file);summary=receipt['summary'];early=receipt['early_measurement']
  market_residual=validate_market(read(collected/case['case_id']/'repetition_01/market_quantity_units.json'))
  rows.append({'index':case['index'],'case_id':case['case_id'],'perturbation':case['perturbation'],'contract_sha256':case['contract_sha256'],'parameters':receipt['parameter_values'],'normalization':summary['final']['normalization'],'asset_price':summary['final']['price'],'reconstructed_market_residual':market_residual,'elapsed_seconds':summary['elapsed_seconds'],'stationary_solves':summary['stationary_solves'],'early_measurement':early,'legacy_stationary_moments':summary['final']['legacy_stationary_moments'],'checkpoint':receipt['checkpoint'],'graph_count':len(receipt['original_graphs']),'numerical_gates_verified':True})
 (root/'raw_case_summary.json').write_text(json.dumps({'status':'complete' if len(rows)==19 else 'partial','cases':rows,'case_count':len(rows),'objective':None,'calibrated_smm':False,'weights_assigned':False,'measurement_note':'Both CPS projections and unavailable housing moments preserved exactly. No candidate selection or SMM.'},indent=2,allow_nan=True)+'\n')
 if rows:
  flat=[]
  for row in rows:
   item={key:row[key] for key in ('index','case_id','contract_sha256','asset_price','reconstructed_market_residual','elapsed_seconds','stationary_solves')}
   item.update({'parameter.'+k:v for k,v in row['parameters'].items()})
   item.update({'normalization.'+k:v for k,v in row['normalization'].items()})
   for projection,observer in row['early_measurement']['fertility'].items():
    item.update({'fertility.'+projection+'.'+k:v for k,v in observer['moments'].items()})
   item.update({'housing_wealth.'+k:v for k,v in row['early_measurement']['housing_wealth']['moments'].items()})
   flat.append(item)
  with (root/'raw_case_summary.csv').open('w',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
 states=data['cases'];errors=[s for s in states if s['status'] in ('failed','collection_validation_failed')]
 receipt={'status':'complete' if len(rows)==19 else 'partial','updated_at_utc':data['checked_at_utc'],'job_id':'17358647','source_commit':'7e872053','all_19_unique_case_file_contracts_verified':True,'case_output_contracts_verified':len(rows),'source_pins_verified':634,'seed_sha256_verified':data['seed_sha256_verified'],'case_count_verified':len(rows),'verified_case_ids':[r['case_id'] for r in rows],'pending_or_running':[s['case_id'] for s in states if s['status']=='pending_or_running'],'failures_or_validation_errors':errors,'current_states':states,'checkpoints_verified':len(rows),'graph_hashes_verified':17*len(rows),'graph_pngs_downloaded':False,'raw_case_summary_sha256':digest(root/'raw_case_summary.json'),'collector_sha256':digest(pathlib.Path(__file__)),'calibrated_smm':False,'weights_assigned':False,'remaining':'none for collection' if len(rows)==19 else 'Await remaining existing jobs; no rerun authorized by collector'}
 if (collected/'pending_scheduling_audit.json').exists():
  receipt['scheduling_audit']=read(collected/'pending_scheduling_audit.json')
 (root/'panel_receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=True)+'\n')
 print(json.dumps({k:receipt[k] for k in ('status','case_count_verified','verified_case_ids','pending_or_running','failures_or_validation_errors')}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--remote',action='store_true');p.add_argument('--skip',default='');p.add_argument('--ingest',type=pathlib.Path);a=p.parse_args()
 if a.remote:remote(set(a.skip.split(','))-{''})
 elif a.ingest:ingest(a.ingest)
 else:raise SystemExit('Specify --remote or --ingest')
