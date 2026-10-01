"""Exact native birth-function toy and decomposition tests; zero household solves."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from forward_flow_decomposition import grouped_flows, three_corner, validate_flows, install_solve_traps

source=Path(__file__).resolve().parents[6]/'code/model/tools/run_e5f_open_population_transition.py'
tree=ast.parse(source.read_text())
function=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='apply_sequential_fertility')
env=dict(np=np,SimpleNamespace=SimpleNamespace,calendar=SimpleNamespace(model=SimpleNamespace(
    independent_child_maturation_active=lambda P:True,get_fecundity_by_age=lambda P:np.array([.5]),
    readiness_settled_state=lambda P:0,birth_destination_child_state=lambda P,m:m+1)))
exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),env)
native=env['apply_sequential_fertility']
P=SimpleNamespace(sequential_births=True,J=1,I=1,n_parity=4,A_f_start=1,A_f_end=1)
shape=(2,1,1,1,1,4,4)
pre=np.zeros(shape);pre[...,0,0]=.25;pre[...,1,0]=.125;pre[...,2,1]=.125
p=np.zeros(shape[:-2]+(4,));p[...,0]=.8;p[...,1]=.2
c=np.zeros(shape[:-2]+(2,2,4));c[...,1,0,:]=.3;c[...,1,1,:]=.4
policy=SimpleNamespace(fert_probs=p,fert2_probs=c)
post,births,_=native(pre,p,P,c)
assert np.isclose(births,.1375)
assert np.isclose(post[...,1,1].sum(),.05) # first births are not exposed to second births
ev=SimpleNamespace(g_pre=pre,g_post_fertility=post,births=births)
saved=dict(first_births=.05,second_births=.0375,third_bin_entries=.05,births=.1375)
validate_flows(ev,saved)
group=grouped_flows(native,pre,policy,P,np.array([-1.,1.]))
assert np.isclose(sum(x['flow'] for x in group.values()),births)
assert np.isclose(sum(x['flow'] for k,x in group.items() if k[3]),births/2)
rows=three_corner([group,{k:dict(flow=v['flow']*2,exposure=v['exposure']) for k,v in group.items()},group])
assert all(np.isclose(r['policy_change']+r['composition_change'],r['total_change']) for r in rows)
assert np.isclose(sum(r['policy_change'] for r in rows),births)
assert np.isclose(sum(r['composition_change'] for r in rows),-births)
try:validate_flows(ev,dict(saved,first_births=.1))
except ValueError:pass
else:raise AssertionError('Saved-flow mismatch accepted')
module=SimpleNamespace(__name__='toy',solve_household=lambda:None,precompute=lambda:True)
assert install_solve_traps(module)==['toy.solve_household'] and module.precompute()
try:module.solve_household()
except RuntimeError:pass
else:raise AssertionError('Solver trap failed')
print(json.dumps(dict(status='passed',household_solves=0,checks=['exact native fecundability and one-birth clock','no newly born child enters same-period second-birth risk','native birth-order replay and mismatch rejection','origin-family/age/wealth groups add to native total','three-corner identity and exposure retention','solver entry trap']),indent=2))
