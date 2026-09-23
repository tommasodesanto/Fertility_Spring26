"""Bounded binding/schema tests; no model equilibrium runs."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT/'code/model/tools'),str(ROOT/'code/model')]
import run_e5f_preference_share_candidate as adapter
from test_e5f_parenthood_utility import old_parameters
from intergen_eqscale_seq_optimized import solver


def parent_module():
    spec=importlib.util.spec_from_file_location('fresh_parent',ROOT/'code/model/tools/e5f_parenthood_utility.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


def test_control_is_exact_existing_adapter():
    p=parent_module(); raw=old_parameters(); first=p.initialize_parenthood_utility(raw)
    second=p.bind_parenthood_utility(first,{'h_P':.8,'beta_annual':.985,'H0':12.})
    assert second.hbar_first_child_jump==.8 and second.beta==.985**4
    assert second.child_room_floor is True and second.delta_alpha==second.delta_alpha_jump==0
    # New adapter does not install any replacement in control mode.
    assert adapter.preference_mode(dict(preference_specification=dict(mapping='floor_control',author_decision='approved_diagnostic'),entry_specification=dict(rule='zero_assets')))=='floor_control'


def test_shares_bind_equivalence_scale_and_native_shared_arrays():
    p=parent_module(); olddomain=p.PARENTHOOD_SEARCH_DOMAIN; adapter.install_share_utility(p)
    assert p.PARENTHOOD_SEARCH_DOMAIN[:8]==olddomain[:8]
    raw=old_parameters(); init=p.initialize_parenthood_utility(raw)
    P=p.bind_parenthood_utility(init,dict(delta_alpha_jump=.10,delta_alpha=.05,beta_annual=.985,H0=12.,kappa_fert=2.))
    assert P.beta==.985**4 and P.rho_hat==P.rho==1/P.beta-1
    assert P.eps_fert==2. and P.eqscale_form=='power' and P.sigma==2.
    assert P.child_room_floor is False
    assert P.c_bar_0==P.c_bar_n==P.hbar_first_child_jump==P.hbar_child_rooms==0
    assert raw.hbar_first_child_jump==.45 # copy discipline
    shared=solver.precompute_shared(P,np.array([0.,1.]))
    alpha=np.asarray(shared.alpha_flat); hb=np.asarray(shared.hb_flat); cb=np.asarray(shared.cb_flat)
    assert np.all(hb==0) and np.all(cb==0)
    assert set(np.round(alpha.reshape(-1),12))=={.733,.583,.533,.483}


@pytest.mark.parametrize('values',[{'h_P':.2},{'delta_alpha':-.01},{'delta_alpha_jump':.251},{'delta_alpha':float('nan')}])
def test_reject_bad_coordinates(values):
    p=parent_module();adapter.install_share_utility(p)
    with pytest.raises(ValueError): p.validate_parenthood_candidate(values)


def test_full_vector_and_stale_floor_fail():
    p=parent_module();adapter.install_share_utility(p)
    with pytest.raises(ValueError):p.validate_parenthood_candidate({},require_complete=True)
    P=p.initialize_parenthood_utility(old_parameters());P.hbar_first_child_jump=.1
    with pytest.raises(ValueError):p.bind_parenthood_utility(P,{})


def test_generated_contracts_preserve_targets_and_numeric_gates(tmp_path):
    source=ROOT/'tmp/earnings_wealth_direct_period_20260922_v5/plan.json'
    if not source.exists(): pytest.skip('Optional retained V5 fixture absent')
    plan=json.loads(source.read_text());plan['preference_specification']={'mapping':'child_dependent_shares','author_decision':'approved_diagnostic'}
    for name,module in [('accounting',adapter.earnings.accounting),('income',adapter.earnings.income),('period_income',adapter.earnings.period_income)]: plan['files'][name]=adapter.pin(module.__file__)
    plan['structural_parameters'].pop('h_P');plan['structural_parameters'].update(delta_alpha=.05,delta_alpha_jump=.1)
    result=adapter.prepare_plan(plan,tmp_path/'generated')
    adapter.verify_plan(result)
    oldrun=adapter.earnings.read(plan['files']['run_contract']['path']);newrun=adapter.earnings.read(result['files']['run_contract']['path'])
    oldobj=adapter.earnings.read(oldrun['working_objective']['path']);newobj=adapter.earnings.read(newrun['working_objective']['path'])
    assert {k:v for k,v in oldobj.items() if k!='parameter_restrictions'}=={k:v for k,v in newobj.items() if k!='parameter_restrictions'}
    assert len(newobj['parameter_restrictions'])==10
    assert result['objective_canonical_sha256']!=plan['objective_canonical_sha256']
    wrapper=Path(result['files']['wrapper']['path']).read_text()
    assert 'len(parameters)==19' in wrapper
    for key in ('scorer','validator'):
        compile(Path(newrun[key]['path']).read_text(),newrun[key]['path'],'exec')
    compile(wrapper,result['files']['wrapper']['path'],'exec')
    before=Path(oldrun['validator']['path']).read_text();after=Path(newrun['validator']['path']).read_text()
    # Every existing gate remains; only fixed-zero assertion extended and coordinate globals appended.
    assert after.startswith(before.replace("{'hbar_child_rooms':0.,'payroll_tax'", "{'hbar_first_child_jump':0.,'hbar_child_rooms':0.,'payroll_tax'"))
    result['target_system_sha256']='0'*64
    with pytest.raises(ValueError,match='Target/weight'):adapter.verify_plan(result)


@pytest.mark.parametrize('mode',['floor_control','child_dependent_shares'])
def test_native_zero_solve_preflight(tmp_path,mode):
    source=ROOT/'tmp/earnings_wealth_direct_period_20260922_v5/plan.json'
    if not source.exists(): pytest.skip('Optional retained V5 fixture absent')
    plan=json.loads(source.read_text())
    for name,module in [('accounting',adapter.earnings.accounting),('income',adapter.earnings.income),('period_income',adapter.earnings.period_income)]:plan['files'][name]=adapter.pin(module.__file__)
    plan['preference_specification']={'mapping':mode,'author_decision':'approved_diagnostic'}
    if mode=='child_dependent_shares':
        plan['structural_parameters'].pop('h_P');plan['structural_parameters'].update(delta_alpha=.05,delta_alpha_jump=.10)
    plan=adapter.prepare_plan(plan,tmp_path/'inputs')
    planpath=tmp_path/'plan.json';adapter.earnings.write(planpath,plan)
    adapter.verify_plan(plan)
    case=next(c for c in plan['cases'] if c['arm']=='literature_income_purchase')
    result=adapter.run_case(planpath,plan,case['arm'],tmp_path/'output',case['repetitions'],True)
    assert result['household_solves']==0 and result['status']=='preflight_passed'
