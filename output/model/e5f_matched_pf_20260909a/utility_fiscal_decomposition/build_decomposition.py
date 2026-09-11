"""Verify the original receipts and tabulate the fixed-psi four-cell comparison."""
import csv,hashlib,json,math,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
DATA=ROOT/'collected'
CASES=('old_old','old_balanced','new_old','new_balanced')
NAMES=('beta_annual','kappa_fert','kappa_fert_continuation','chi','H0','theta0','theta1','first_birth_fixed_cost','h_P')
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def write(p,value): p.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def csvout(path,rows):
 with path.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def finite(value):
 return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)
manifest=read(DATA/'collection_manifest.json'); contract=read(ROOT/'contract.json'); submission=read(ROOT/'submission_receipt.json')
assert digest(ROOT/'contract.json')==submission['contract_sha256']
for name,pin in manifest['files'].items():
 assert digest(DATA/name)==pin['sha256'],name
summaries={}; params={}; metrics={}; case_rows=[]; parameter_rows=[]
for case in CASES:
 path=DATA/case; summary=read(path/'summary.json'); final=summary['final']; summaries[case]=summary
 assert summary['status']=='passed_initial_candidate_loop' and summary['normalized'] is False
 assert summary['calibrated_smm'] is False and summary['perfect_foresight_solved'] is False
 assert summary['stationary_solves']==summary['repetitions']==1
 assert final['checkpoint_sha256']==manifest['checkpoints'][case]['sha256']
 assert read(path/'contract.json')['contract_sha256']==submission['contract_sha256']
 seed=read(path/'seed_mapping.json')
 assert seed['initial_psi']==contract['frozen_psi_child']
 assert seed['source_files_verified']==625
 with (path/'repetition_01/parameters.csv').open() as f: rows=list(csv.DictReader(f))
 params[case]={row['parameter']:float(row['estimate']) for row in rows}
 assert params[case]['psi_child']==contract['frozen_psi_child']
 assert params[case]['payroll_tax']==0.179 and params[case]['housing_supply_elasticity']==0.63
 if case.startswith('new_'): assert params[case]['hbar_child_rooms']==0.0
 for row in rows: parameter_rows.append({'case':case,**row})
 accounts=final['fiscal']['actual_accounts']
 if case.endswith('_balanced'):
  assert final['fiscal']['fiscal_gate'] and abs(accounts['scaled_pension_budget_residual'])<=1e-6
 values={'asset_price':final['price']}
 values.update({'legacy_moment.'+k:v for k,v in final['legacy_stationary_moments'].items()})
 values.update({'fiscal.'+k:v for k,v in accounts.items()})
 quantity=read(path/'repetition_01/market_quantity_units.json')
 values['quantity.adult_households']=quantity['adult_households']
 values.update({'quantity.'+k:v for k,v in quantity['aggregate_quantities'].items()})
 metrics[case]=values
 ge=read(path/'repetition_01/stationary_solves.json')
 assert len(ge)==1 and ge[0]['status']=='completed'
 case_rows.append({'case':case,'utility':'parenthood_only' if case.startswith('new_') else 'old_jump_plus_slope',
  'pension_rule':'actual_budget_balanced' if case.endswith('_balanced') else 'old_unbalanced_diagnostic',
  'job_id':submission['job_id']+'_'+str(CASES.index(case)),'elapsed_seconds':summary['elapsed_seconds'],
  'GE_seconds':ge[0]['seconds'],'asset_price':final['price'],'psi_child':params[case]['psi_child'],
  'pension_period':accounts['pension_period_units'],'fiscal_scaled_residual':accounts['scaled_pension_budget_residual'],
  'GE_market_error':ge[0]['market_error'],'checkpoint_path':manifest['checkpoints'][case]['path'],
  'checkpoint_sha256':final['checkpoint_sha256'],'original_graph_count':len(manifest['original_graphs'][case]),
  'calibrated_smm':False,'observer_status':'legacy stationary observers; not certified early targets'})
for name in NAMES:
 assert len({params[c][name] for c in CASES})==1,(name,[params[c][name] for c in CASES])
assert params['old_old']['hbar_child_rooms']==params['old_balanced']['hbar_child_rooms']
assert params['old_old']['hbar_child_rooms']>0
rows=[]; unavailable=[]
for name in sorted(set.union(*(set(v) for v in metrics.values()))):
 values={case:metrics[case].get(name) for case in CASES}
 if not all(finite(v) for v in values.values()):
  unavailable.append({'metric':name,'reason':'not a finite scalar in all four cases'})
  continue
 a,b,c,d=(float(values[case]) for case in CASES)
 rows.append({'metric':name,'observer_status':'legacy stationary diagnostic' if name.startswith('legacy_moment.') else 'raw model accounting diagnostic',
  **values,'utility_effect_old_pension':c-a,'utility_effect_balanced_pension':d-b,
  'pension_effect_old_utility':b-a,'pension_effect_new_utility':d-c,'interaction':d-b-c+a,
  'combined_change':d-a})
 for key in ('utility_effect_old_pension','utility_effect_balanced_pension','pension_effect_old_utility','pension_effect_new_utility','interaction','combined_change'):
  assert math.isfinite(rows[-1][key])
 assert math.isclose((c-a)+(d-c),d-a,rel_tol=1e-12,abs_tol=1e-12)
write(ROOT/'raw_decomposition.json',{'classification':'frozen-preference stationary diagnostic; not calibrated SMM or perfect foresight',
 'source_commit':'c6dd3508','job_id':submission['job_id'],'contract_sha256':submission['contract_sha256'],
 'source_artifact_count_verified':len(manifest['files']),'unchanged_structural_coordinates':{name:params['old_old'][name] for name in NAMES},
 'frozen_psi_child':contract['frozen_psi_child'],'case_results':case_rows,'raw_comparisons':rows,'unavailable_or_nonscalar_metrics':unavailable,
 'difference_definitions':{'utility_effect_old_pension':'new_old - old_old','utility_effect_balanced_pension':'new_balanced - old_balanced',
 'pension_effect_old_utility':'old_balanced - old_old','pension_effect_new_utility':'new_balanced - new_old','interaction':'new_balanced - old_balanced - new_old + old_old'},
 'measurement_warning':'No early target values, gaps, weights or objective are applied. Legacy stationary fertility timing and parent-group observers need the separate early-contract review. Old-pension cases are intentionally fiscally invalid and cannot be selected benchmarks.'})
csvout(ROOT/'raw_decomposition.csv',rows);csvout(ROOT/'case_summary.csv',case_rows);csvout(ROOT/'parameters_all_cases.csv',parameter_rows)
write(ROOT/'original_graphs_receipt.json',manifest['original_graphs'])
print(json.dumps({'cases':4,'raw_scalar_metrics':len(rows),'original_graphs':68,'light_artifacts_verified':len(manifest['files']),'summaries':case_rows},indent=2))
