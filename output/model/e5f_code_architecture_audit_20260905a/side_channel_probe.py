from pathlib import Path
import sys,json,hashlib
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
root=Path.cwd(); frozen=root/'tmp/e5f_overnight_20260905a/frozen_production/code/model';sys.path[:0]=[str(frozen/'tools'),str(frozen)]
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
from intergen_eqscale_seq_optimized import solver
calendar.model=solver;calendar.apply_fertility=transition.apply_sequential_fertility
P=SimpleNamespace(J=1,I=1,n_house=0,n_parity=4,n_child_states=4,A_f_start=1,A_f_end=1,sequential_births=True,child_state_mode='independent_count',readiness_gate_enabled=False,fecundity_omega1=0.,fecundity_omega2=0.,age_start=18,da=4.,user_cost_rate=.1,H0=np.array([1.]),r_bar=.1,xi_supply=np.array([1.]))
g=np.zeros((2,1,1,1,1,4,4));g[0,0,0,0,0,1,1]=1.
fp=np.zeros((2,1,1,1,1,4));fp[...,0]=1.
pol=calendar.PolicyBundle(V=np.ones_like(g),c_pol=np.ones_like(g),hR_pol=np.ones_like(g),bp_pol=np.zeros_like(g),tenure_choice=np.zeros_like(g,dtype=int),tenure_probs=None,loc_probs=np.ones((2,1,1,1,1,1,4,4)),fert_probs=fp,fert_value=np.zeros((2,1,1,1,1)),price=np.array([1.]),maps=SimpleNamespace(lmm_idx=None,lmm_wt=None,tmx_idx=None,tmx_wt=None))
cont=np.zeros((2,1,1,1,1,2,2,4));P._fert2_probs=cont.copy()
with patch.object(calendar,'gate_pre_fertility_distribution',side_effect=lambda g,*args:(g.copy(),0.)),patch.object(solver,'realize_current_cross_section',side_effect=lambda g,*args,**kwargs:g.copy()):
 a=calendar.evaluate_period(np.array([1.]),g,P,np.array([0.,1.]),SimpleNamespace(),calendar.SolveCounter(),supplied_policy=pol)
 P._fert2_probs[...,1,0,1]=1.
 b=calendar.evaluate_period(np.array([1.]),g,P,np.array([0.,1.]),SimpleNamespace(),calendar.SolveCounter(),supplied_policy=pol)
assert a.births==0. and b.births==1.
receipt={'test':'identical supplied PolicyBundle depends on unrelated last-solve continuation side channel','source':'frozen production','source_files':{str(p.relative_to(frozen)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(calendar.__file__),Path(transition.__file__),Path(solver.__file__)]},'first_evaluation_births':a.births,'second_evaluation_births':b.births,'same_policy_object':a.policy is b.policy,'same_inherited_distribution':np.array_equal(a.g_pre,b.g_pre),'same_price':np.array_equal(a.policy.price,b.policy.price),'changed_only_parameter_field':'P._fert2_probs','mocked_operations':['dead-mass gate identity','current-tenure realization identity'],'production_effect_demonstrated':False,'interpretation':'Confirmed non-self-contained policy API. Actual fertility operator and evaluate_period executed; no economic solve or production result changed.'}
(root/'output/model/e5f_code_architecture_audit_20260905a/side_channel_probe.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
