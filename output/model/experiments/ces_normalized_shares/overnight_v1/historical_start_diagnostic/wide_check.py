import json,time,numpy as np
from pathlib import Path
from experiments.ces_normalized_shares import adapter
from production import equilibrium
root=Path('/work/deployment'); packet=json.loads((root/'followup_tools/historical_start_point.json').read_text())
out=Path('/work/results/run');out.mkdir(parents=True,exist_ok=False)
end=time.time()+1500
start=6.318956013969669+(10.110329622351472-6.318956013969669)*0.1025104211342307/(0.1025104211342307+0.16786110455107384)
(out/'plan.json').write_text(json.dumps(dict(packet,diagnostic_only=True,no_search=True,numerical_price_caps_origin_scale=2,price_start=start,deadline_epoch=end,maximum_lifecycle_solves=32,unchanged='All economic inputs, preference normalization and scoring unchanged; only numerical bracketing origin and initial trial price differ.'),indent=2)+'\n')
def wider_solver(*args,**kwargs):
 original=equilibrium.build_context
 def build(*a,**k):
  c=original(*a,**k)
  c['q_ref']=2*np.asarray(c['q_ref']).copy()
  return c
 equilibrium.build_context=build
 kwargs['price_start']=start
 try:return equilibrium.solve_stationary_ge(*args,**kwargs)
 finally:equilibrium.build_context=original
P,grid=adapter.load_inputs(packet['parameters'])
evaluate=adapter.make_evaluator(out,'historical_start_wider_caps',P,grid,end,solver=wider_solver,target_fingerprint=packet['target_fingerprint'],weight_fingerprint=packet['weight_fingerprint'])
result=evaluate('historical_start',packet['parameters'],end)
assert result['effective_input_fingerprint']=='316f9b66ac13bd9f77110805d48723141249668a1e9949afa6d3f32062dc571e'
(out/'completed.json').write_text(json.dumps(dict(result,diagnostic_only=True,no_search=True,numerical_price_caps_origin_scale=2,source=packet['source']),indent=2)+'\n')
print(json.dumps({k:result.get(k) for k in ('status','reason','loss','price','H0_derived','lifecycle_solves','report')}))
