#!/usr/bin/env python3
"""Build the lead-specified diagnostic plan from target-architecture identity."""
import csv,hashlib,importlib.util,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5];OUT=Path(__file__).resolve().parent
CODE=ROOT/'code/model/experiments/transition_readiness'
def pin(p):return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def load(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
controller=load('current_floor_controller',CODE/'one_shock_floor.py')
pre=json.loads((OUT/'native_import_preflight.json').read_text())
assert pre['native_setup_verified'] and pre['policy_calls']==0
handoff=ROOT/'output/model/transition_readiness_v1/current_floor_handoff/handoff.json'
h=json.loads(handoff.read_text())
with (handoff.parent/'parameters.csv').open(newline='') as f:
 psi=next(float(r['estimate']) for r in csv.DictReader(f) if r['parameter']=='psi_child')
contract=json.loads((OUT.parent/'controller/target_contract.json').read_text())
plan=dict(schema='current_floor_one_permanent_v1',kind='one_permanent',start_year=2007,mode='diagnostic',horizons=[6,8],
 fiscal_relaxation_authorized=True,gates=controller.GATES.copy(),identity=pre['identity'],handoff=pin(handoff),
 source_files=dict(controller=pin(CODE/'one_shock_floor.py'),runtime=pin(CODE/'floor_runtime.py'),
 original_estimator=pin(CODE/'pinned_tools/run_e5f_preference_estimation.py'),original_fitter=pin(CODE/'pinned_tools/e5f_preference_shock_fit.py')),
 target_contract=contract,execution_enabled=True,
 budget=dict(total_seconds=5400,seed_seconds=3600,candidate_seconds=3000,endpoint_seconds=1800,mapping_seconds=900,path_seconds=2400,render_seconds=300,maximum_policy_calls=260),
 seed=dict(horizon=12,perturbed_date=5,log_step=1e-5),initial_psi=psi,native_smoke_psi=psi*1.001,psi_bound_ratios=[.01,2.],
 standard_plot_names=sorted(h['selected_repeat_verification']['standard_plot_sha256']),
 fit=dict(max_evaluations=9,log_difference_step=.01,fertility_tolerance=.005,max_log_step=.15,damping=.7,max_condition_number=1e8,worsening_factor=1.5,reproduction_tolerance=1e-8),
 path=dict(max_evaluations=6,price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],max_log_step=.15,damping=.7),
 endpoint=dict(max_evaluations=18,price_bound_ratios=[.05,20.],max_log_step=.15,damping=.7,slope=1.),
 economic_disclosure=dict(status='author-selected current-floor transition experiment',reference='verified chain0/case0006_nm',changes_relative_to_block0506=['parenthood housing floor','current recalibrated preferences and dispersion','empirical nonnegative mean-preserving entrant wealth','corrected zero unsecured credit','current static-elastic housing supply and physical population scaling'],stationary_price=h['selected_point']['price'],stationary_population=h['selected_point']['population'],dated_fiscal_tolerance=2e-5,stationary_fiscal_tolerance=1e-6,full_104_128_not_yet_verified=True,no_automatic_fallback=True))
receipt=controller.preflight(plan)
(OUT/'smoke_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
(OUT/'smoke_plan_preflight.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
