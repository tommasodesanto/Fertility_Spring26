import ast,csv,hashlib,json,pathlib,shutil
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent
RUN=HERE.parent; C=HERE/'candidate_0001'
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(pathlib.Path(p).read_text())
planpath=RUN/'executed/fit_plan.json';plan=read(planpath);inventory=read(RUN/'executed/inventory.json')
complete=read(C/'complete.json');checkpoint=read(C/'state_2023_checkpoint/checkpoint_receipt.json')
assert checkpoint['identity']==plan['identity']
assert sha(C/'state_2023_checkpoint/actual_2023.pkl.gz')==checkpoint['export']['sha256']
assert checkpoint['export']==read(C/'state_2023_checkpoint/actual_2023.json')
assert checkpoint['state_experiment_ready'] is False and complete['shock_fit_complete'] is False
assert checkpoint['export']['exact_native_state'] and checkpoint['export']['no_rescaling'] and checkpoint['export']['both_native_queues']
shutil.copy2(planpath,HERE/'executed_fit_plan.json')
source=RUN/'executed/source/code/model/experiments/transition_readiness/pinned_tools/e5f_preference_shock_fit.py'
assert sha(source)==inventory['files']['code/model/experiments/transition_readiness/pinned_tools/e5f_preference_shock_fit.py']
for key in ('blocks','annual'):
 assert sha(plan['target_contract'][key]['path'])==plan['target_contract'][key]['sha256']
# Execute only the archived deterministic table function, not its model-facing imports.
function=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='fit_rows')
ns={'np':np};exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),ns)
records=[]
for h in (12,16):
 directory=C/f'horizon_{h:03d}';root=read(directory/'root.json');m1=read(directory/'map_001/native_record.json');m2=read(directory/'map_002/native_record.json');audits=read(directory/'map_002/dated_audits.json')
 assert root['converged'] and all(root['gates'].values()) and all(m2['gates'].values()) and len(audits)==h
 assert all(all(a['gates'].values()) for a in audits)
 assert root['final_reproduction_max_abs']==0.
 gap=max(abs(x-y) for field in ('market_residual','fiscal_residual') for x,y in zip(m1[field],m2[field]))
 assert gap==0.
 records.append(dict(horizon=h,root_status=root['status'],root_evaluations=root['evaluations'],root_gates=root['gates'],
   native_gates=m2['gates'],dated_accounting_rows=len(audits),all_dated_accounting_pass=True,
   maximum_market_residual=max(map(abs,m2['market_residual'])),maximum_fiscal_residual=max(map(abs,m2['fiscal_residual'])),
   exact_replay_gap=gap,mass_error=m2['mass_error'],policy_reproduction_error=m2['backward_forward_policy_error'],
   projection_mass=m2['projection_mass'],native_policy_calls_map1=m1['policy_calls'],native_policy_calls_map2=m2['policy_calls']))
models=[r['period_tfr_topcode_adjusted'] for r in m2['fertility'][:4]]
assert models==complete['payload']['models']
rows=ns['fit_rows'](plan['target_contract']['rows'],models,'one_permanent')
assert sum(r['loss_contribution'] for r in rows)==complete['loss_contribution']
(HERE/'target_fit.json').write_text(json.dumps(rows,indent=2)+'\n')
columns=['moment','decision_year','birth_year_start','birth_year_end','target','model','gap','weight','loss_contribution','role']
with (HERE/'target_fit.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
psi=complete['psi'];lo,hi=[plan['initial_psi']*x for x in plan['psi_bound_ratios']]
parameter=dict(parameter='permanent_2007_psi_child',estimate=psi,lower=lo,upper=hi,near_bound=min(psi-lo,hi-psi)<=.01*(hi-lo),
 status='baseline candidate only; fitted estimate not yet obtained',reference_estimate=plan['initial_psi'])
with (HERE/'shock_parameter.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(parameter));writer.writeheader();writer.writerow(parameter)
handoffpath=pathlib.Path(plan['handoff']['path']);assert sha(handoffpath)==plan['handoff']['sha256']
handoff=read(handoffpath);table=pathlib.Path(handoff['tables']['parameters_csv']['path']) if 'path' in handoff['tables']['parameters_csv'] else pathlib.Path('output/model/transition_readiness_v1/current_floor_handoff/parameters.csv')
table=table if table.is_absolute() else handoffpath.parent/table
assert sha(table)==handoff['tables']['parameters_csv']['sha256'];shutil.copy2(table,HERE/'baseline_parameters_all31.csv')
endpoint=read(C/'horizon_012/endpoint/latest_completed.json')
assert endpoint['accounting_valid'] and abs(endpoint['renewal_residual'])<=plan['gates']['stationary_renewal_tolerance']
summary=dict(endpoint=endpoint,classification='Unchanged-preference baseline candidate; not a fitted transition or experiment-ready state',candidate=1,
 mode=plan['mode'],shock_fit_complete=False,production_ready=False,state_experiment_ready=False,full_path_certified=False,
 psi=psi,weighted_loss=complete['loss_contribution'],target_fit_rows=4,parameter=parameter,horizons=records,
 state_physical_horizon_stable=checkpoint['state_physical_horizon_stable'],state_value_integrity=checkpoint['state_value_integrity'],
 terminal_passes=checkpoint['terminal_passes'],native_state=checkpoint['export'],
 local_native_state=dict(path=str(C/'state_2023_checkpoint/actual_2023.pkl.gz'),sha256=sha(C/'state_2023_checkpoint/actual_2023.pkl.gz'),bytes=(C/'state_2023_checkpoint/actual_2023.pkl.gz').stat().st_size),
 evidence_pins=dict(executed_fit_plan=sha(planpath),archived_fit_rows_source=sha(source),checkpoint_receipt=sha(C/'state_2023_checkpoint/checkpoint_receipt.json'),
    handoff=sha(handoffpath),baseline_parameter_table=sha(table)),
 unresolved=['Shock fitting remains in progress; this candidate keeps baseline preferences.','Diagnostic horizons12/16 do not satisfy production104/128 certification.','Checkpoint receipt explicitly keeps state_experiment_ready=false.','Visual diagnostic packet has not been collected or reviewed in this evidence-only pass.','Model household-rate analogue differs from published female-exposure TFR.'],
 native_model_calls_during_collection=0)
(HERE/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:summary[k] for k in ('weighted_loss','horizons','local_native_state')},indent=2))
