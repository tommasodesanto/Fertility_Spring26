from pathlib import Path
import sys,json
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
r=Path.cwd(); f=r/'tmp/e5f_overnight_20260905a/frozen_production/code/model';sys.path[:0]=[str(f/'tools'),str(f)]
import run_dynamic_population_transition as c
calls=[]
def fake(price,*args,**kwargs):
 calls.append(float(price[0]));return SimpleNamespace(demand_by_loc=np.array([2.]),supply_by_loc=np.array([1.]),relative_market_residual=1.)
with patch.object(c,'evaluate_period',side_effect=fake):
 try:c.clear_scalar_housing_market(np.ones(1),1.,SimpleNamespace(I=1,p_min=.1,p_max=1.1),np.array([0.,1.]),SimpleNamespace(),c.SolveCounter(),1e-4,2,c.HousingSupplyRule('fixed-stock',1.,1.,0.))
 except RuntimeError as exc:error=str(exc)
assert len(calls)==19 and calls.count(1.1)==18
out={'source':'frozen production','test':'unbracketed positive excess at upper price bound','evaluate_period_calls':len(calls),'distinct_prices':len(set(calls)),'same_upper_bound_evaluations':calls.count(1.1),'prices':calls,'raised':error,'model_solves_performed':0,'interpretation':'Mocked residual demonstrates redundant boundary evaluations; actual runtime occurrence not established.'}
(r/'output/model/e5f_code_architecture_audit_20260905a/bracket_probe.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
