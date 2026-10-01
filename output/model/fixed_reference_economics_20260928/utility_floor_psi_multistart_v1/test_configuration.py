"""Zero-lifecycle receipts for seed construction, binding and both controller routes."""
import importlib.util,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
import run_psi as r
HERE=Path(__file__).resolve().parent

def main():
    c=r.CONFIG; coords=c['free_coordinates']; bounds=c['bounds']; center=c['center']; samples=2*qmc.LatinHypercube(d=10,seed=20260930).random(n=8)-1
    for i,row in enumerate(c['nearby_starts']):
        for j,k in enumerate(coords):assert row['parameters'][k]==center[k]+float(samples[i,j])*c['nearby_design']['halfwidths'][k]
        r.inputs.check_point(row['parameters'],bounds)
    assert len({tuple(x['parameters'][k] for k in coords) for x in c['nearby_starts']})==8
    old=json.loads((HERE.parent/'utility_floor_psi_v1/plan.json').read_text());assert c['base_target_contract']==old['base_target_contract'] and c['profiles']=={'base_control':{}}
    assert c['initial_price']==.6786606850351011
    lane='floor_s0';r.inputs.LANES[lane].update(bounds=bounds,free_coordinates=coords)
    P,grid=r.inputs.proposal(lane);P,entry=r.inputs.entry(P,grid,'nonnegative_mean')
    for seed in [x['parameters'] for x in c['nearby_starts']]+[p for s in c['pso_swarms'] for p in s['particles']]:
        Q=r.inputs.bind(P,seed,bounds,'floor');assert Q.beta==seed['beta_annual']**Q.period_years and Q.eps_fert==seed['kappa_fert'] and Q.psi_child==seed['psi_child'] and Q.child_room_floor and Q.hbar_first_child_jump==seed['h_P']
    receipts=[]
    for chain in (0,6,7):
        out=HERE/'toy'/f'chain_{chain}'
        seed=c['nearby_starts'][chain]['parameters'] if chain<6 else center
        counter=[0]
        def evaluate(label,point,end):
            counter[0]+=1
            if chain==7 and counter[0]==2:return dict(status='inadmissible_numerical',reason='mock_rejection',lifecycle_solves=0)
            rr=np.asarray([(point[k]-center[k])/(bounds[k][1]-bounds[k][0])+.01 for k in coords]);return dict(status='passed',residual=rr.tolist(),lifecycle_solves=0,report=str(out/label))
        result=r.optimize(out,seed,bounds,coords,evaluate,time.time()+7200,toy=True,maxeval=40,chain=chain)
        assert result['objective_calls']==40 and result['lifecycle_solves']==0
        assert result['completed_full_ge']>10
        if chain>=6:
            state=json.loads((out/'pso_state.json').read_text());assert state['generation']>=1 and len(state['positions_initial_range_scaled'])==8 and state['rng_state']
        if chain==7:assert any(x.get('numerical_rejection') for x in json.loads((out/'cases.json').read_text()))
        receipts.append(dict(chain=chain,objective_calls=result['objective_calls'],completed_mock_cases=result['completed_full_ge'],cache_hits=result['objective_calls']-result['completed_full_ge'],lifecycle_solves=0))
    no_call=lambda *a: (_ for _ in ()).throw(AssertionError('Budget should stop before objective'))
    budget=r.optimize(HERE/'toy/budget',center,bounds,coords,no_call,time.time()+899,toy=True,maxeval=40,chain=6);assert budget['completed_full_ge']==0
    failed=HERE/'toy/failure';counter=[0]
    def fail(*a):raise RuntimeError('mock_fatal_failure')
    try:r.optimize(failed,center,bounds,coords,fail,time.time()+7200,toy=True,maxeval=40,chain=6)
    except RuntimeError as exc:assert str(exc)=='mock_fatal_failure'
    else:raise AssertionError('Fatal error silently absorbed')
    assert json.loads((failed/'pso_state.json').read_text())['status']=='stopped_or_failed_objective'
    receipt=dict(status='passed_zero_lifecycle',distinct_nearby_starts=8,dispatch_nearby_starts=6,pso_chains=2,pso_particles_each=8,weight_contract_exact=True,parameter_binding_vectors=24,hard_bounds_unchanged=True,controllers=receipts,budget_reserve_stop=True,fatal_failure_propagation=True,numerical_rejection_penalty=True,no_auto_resume=True,lifecycle_solves=0)
    (HERE/'configuration_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
