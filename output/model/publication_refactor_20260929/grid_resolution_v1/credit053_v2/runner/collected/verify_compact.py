"""Validate collected hashes, complete fits/parameters, plots and recorded gates."""
from pathlib import Path
import csv,json,hashlib,math
HERE=Path(__file__).resolve().parent;root=HERE/'full'
receipt=json.loads((HERE/'remote_hash_receipt.json').read_text())
for rel,row in receipt['files'].items():assert hashlib.sha256((HERE/rel).read_bytes()).hexdigest()==row['sha256']
assert all(receipt['source_inventory_verified'].values()) and len(receipt['source_inventory_verified'])==81
assert receipt['archive_sha256']=='ea7a8deb347d54a5c57bd4e4580ccbf23c0cbaefc801e5416408c929d960813c'
parameters=[];arms={}
for arm in ['control_160x15','proposal_120x9']:
    prefix=root/arm;completed=json.loads((prefix/'completed.json').read_text());assert completed['result']['status']=='passed';assert completed['lifecycle_solves']==7
    tables=[]
    for final in ['selected_root','selected_repeat_final']:
        folder=prefix/'phase_b_ge'/final
        with (folder/'target_fit.csv').open(newline='') as f:fits=list(csv.DictReader(f))
        with (folder/'parameters.csv').open(newline='') as f:params=list(csv.DictReader(f))
        assert len(fits)==14 and len(params)==31;tables.append((fits,params))
        for r in fits:
            if r['target'] and r['model']:assert math.isclose(float(r['gap']),float(r['model'])-float(r['target']),abs_tol=1e-12)
            if r['weight']:assert math.isclose(float(r['loss_contribution']),float(r['weight'])*float(r['gap'])**2,rel_tol=2e-12,abs_tol=1e-12)
        gates=json.loads((folder/'gates.json').read_text());assert gates['feasibility_projection_mass']==0 and gates['household_budget']['budget_excess_mass']<=2e-10
        purchase=gates['purchase'];assert purchase['maximum_occupied_transaction_wealth_error']<=1e-9
        for key,value in purchase.items():
            if key.endswith('violation_mass') or key in ['transaction_outside_grid_mass','negative_estate_exposure_mass','saving_outside_grid_mass']:assert value<=2e-10
        assert gates['estate']['status']=='funded';assert gates['policy_arrays']['occupied_negative_steps']==0
        for p in gates['policy_arrays']['probabilities'].values():assert not p['nonfinite'] and 0<=p['minimum']<=p['maximum']<=1
        assert abs(gates['fiscal_certificate']['actual_accounts']['scaled_pension_budget_residual'])<=1e-6
        for key in ['stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1','mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error']:assert abs(gates['stationary_operator'][key])<=5e-9
        assert abs(gates['stationary_operator']['zero_entry_mass_accounting_residual'])<=2e-8
        closure=json.loads((folder/'closure.json').read_text());assert abs(closure['renewal_residual'])<=1e-6 and abs(closure['absolute_housing_residual'])<=1e-10 and closure['d_bar']==.53
    assert tables[0]==tables[1];parameters.append(tables[0][1])
    selected={Path(rel).name:row['sha256'] for rel,row in receipt['files'].items() if '/'+arm+'/phase_b_ge/selected_root/standard_diagnostics/' in rel}
    repeated={Path(rel).name:d for rel,d in receipt['repeat_plot_hashes'].items() if '/'+arm+'/' in rel}
    assert len(selected)==len(repeated)==17 and selected==repeated
    time=json.loads((prefix/'external_workflow_timing.json').read_text())['seconds']
    arms[arm]=dict(lifecycle_solves=7,external_workflow_seconds=time,scalar_loss=sum(float(r['loss_contribution'] or 0) for r in tables[0][0]),closure=completed['result']['selected'])
differences=[(a['parameter'],a['estimate'],b['estimate']) for a,b in zip(*parameters) if a!=b]
assert {row[0] for row in differences}=={'wealth_grid_nodes','income_states'}
speed=arms['control_160x15']['external_workflow_seconds']/arms['proposal_120x9']['external_workflow_seconds']
result=dict(status='completed_recorded_scientific_gates_and_compact_receipts_verified',arms=arms,source_files_verified=81,downloaded_files_verified=len(receipt['files']),selected_png_files=34,repeat_png_hashes=34,exact_repeat_fit_parameter_tables_and_plot_hashes=True,native_array_repeat_check='Pinned driver enforces exact native solution/shared array equality before completed status; large arrays not independently downloaded',different_parameter_rows=differences,speed_ratio=speed,workflow_time_reduction=1-1/speed,full_behavior_equivalence_certified=False,policy_occupancy_evidence='policy_extrema_occupancy.json; saved beginning distribution, not postdecision sol.g',additional_model_solves=0)
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(speed_ratio=speed,time_reduction=result['workflow_time_reduction'],losses={k:v['scalar_loss'] for k,v in arms.items()},verified_files=len(receipt['files']))))
