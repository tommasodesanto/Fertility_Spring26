import ast,csv,hashlib,json,pathlib,shutil
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent;C=HERE/'candidate_0004';RUN=HERE.parent
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(pathlib.Path(p).read_text())
planpath=pathlib.Path('output/model/transition_readiness_v1/current_floor_preparation/deployment/fit_v2_plan.json')
assert sha(planpath)=='6b6b28afbf798562f2f28fd8c38eb33d844a36363bc53c45d4b35b87aeb34945'
plan=read(planpath);shutil.copy2(planpath,HERE/'executed_fit_plan.json')
source=RUN/'executed/source/code/model/experiments/transition_readiness/pinned_tools/e5f_preference_shock_fit.py'
assert sha(source)==read(RUN/'executed/inventory.json')['files']['code/model/experiments/transition_readiness/pinned_tools/e5f_preference_shock_fit.py']
for key in ('blocks','annual'):assert sha(plan['target_contract'][key]['path'])==plan['target_contract'][key]['sha256']
function=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='fit_rows')
ns={'np':np};exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),ns)
checkpoint=read(C/'state_2023_checkpoint/checkpoint_receipt.json')
assert checkpoint['identity']==plan['identity']
assert sha(C/'state_2023_checkpoint/actual_2023.pkl.gz')==checkpoint['export']['sha256']
assert sha(C/'state_2023_checkpoint/checkpoint_receipt.json')=='a3e996e74c0de2d53b1e66e7f9de4901cb2cb6950a99ed2a9d723a9539db1de1'
assert checkpoint['shock_fit_complete'] is False and checkpoint['state_experiment_ready'] is False and checkpoint['full_path_certified'] is False
assert not (C/'complete.json').exists()
comparison=read(C/'horizon_comparison.json');state=read(C/'state_2023_horizon_comparison.json')
assert comparison['passed'] is False and state['passed'] is False
summary=[];columns=['moment','decision_year','birth_year_start','birth_year_end','target','model','gap','weight','loss_contribution','role']
for h in (12,16):
 d=C/f'horizon_{h:03d}';root=read(d/'root.json');maps=sorted(d.glob('map_*'));m=read(maps[-1]/'native_record.json');prev=read(maps[-2]/'native_record.json');a=read(maps[-1]/'dated_audits.json')
 assert root['converged'] and all(root['gates'].values()) and all(m['gates'].values())
 assert all(all(x['gates'].values()) for x in a) and len(a)==h
 replay=max(abs(x-y) for key in ('market_residual','fiscal_residual') for x,y in zip(m[key],prev[key]));assert replay<=plan['gates']['final_reproduction_tolerance']
 rows=ns['fit_rows'](plan['target_contract']['rows'],[r['period_tfr_topcode_adjusted'] for r in m['fertility'][:4]],'one_permanent')
 with (HERE/f'target_fit_horizon_{h:03d}.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)
 (HERE/f'target_fit_horizon_{h:03d}.json').write_text(json.dumps(rows,indent=2)+'\n')
 summary.append(dict(horizon=h,models=[r['model'] for r in rows],weighted_loss=sum(r['loss_contribution'] for r in rows),
  root_converged=True,root_evaluations=root['evaluations'],root_gates=root['gates'],mapping_gates=m['gates'],dated_accounting_count=len(a),dated_accounting_pass=True,
  market_maximum=max(map(abs,m['market_residual'])),fiscal_maximum=max(map(abs,m['fiscal_residual'])),replay_gap=replay,
  policy_error=m['backward_forward_policy_error'],projection_mass=m['projection_mass'],mass_error=m['mass_error']))
psi=checkpoint['psi'];lo,hi=[plan['initial_psi']*x for x in plan['psi_bound_ratios']]
param=dict(parameter='permanent_2007_psi_child',estimate=psi,lower=lo,upper=hi,near_bound=min(psi-lo,hi-psi)<=.01*(hi-lo),status='rejected horizon-comparison proposal; shock unfitted',reference_estimate=plan['initial_psi'])
with (HERE/'shock_parameter.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(param));w.writeheader();w.writerow(param)
handoffpath=pathlib.Path(plan['handoff']['path']);assert sha(handoffpath)==plan['handoff']['sha256']
handoff=read(handoffpath);table=pathlib.Path(handoff['tables']['parameters_csv']['path']);table=table if table.is_absolute() else handoffpath.parent/table
assert sha(table)==handoff['tables']['parameters_csv']['sha256'];shutil.copy2(table,HERE/'baseline_parameters_all31.csv')
inventory=read(HERE/'remote_inventory.json')
v=dict(status='rejected_horizon_measurement_comparison',candidate=4,psi=psi,horizons=summary,
 shock_fit_complete=False,state_experiment_ready=False,full_path_certified=False,production_ready=False,
 measurement_horizon_comparison=comparison,state_horizon_comparison=state,terminal_passes=checkpoint['terminal_passes'],
 checkpoint_receipt_sha256=sha(C/'state_2023_checkpoint/checkpoint_receipt.json'),checkpoint_original_flags={k:checkpoint[k] for k in ('shock_fit_complete','one_shock_measurement_certified','horizon_stable','state_experiment_ready','full_path_certified') if k in checkpoint},
 actual_state=dict(path=str(C/'state_2023_checkpoint/actual_2023.pkl.gz'),sha256=sha(C/'state_2023_checkpoint/actual_2023.pkl.gz'),bytes=(C/'state_2023_checkpoint/actual_2023.pkl.gz').stat().st_size),
 files=len(inventory),input_plan_sha256=sha(planpath),archived_fit_rows_sha256=sha(source),
 copied_remote_hashes_verified=True,native_model_calls_during_collection=0,
 unresolved=['Shock not fitted; candidate measurement rejected by unchanged horizon criterion.','Both diagnostic terminal checks fail.','Physical2023 distribution comparison fails.','No production104/128 horizon certification.'],
 table_sha256={p.name:sha(p) for p in HERE.glob('*.csv')})
(HERE/'verification.json').write_text(json.dumps(v,indent=2)+'\n')
print(json.dumps({k:v[k] for k in ('files','psi','horizons','actual_state')},indent=2))
