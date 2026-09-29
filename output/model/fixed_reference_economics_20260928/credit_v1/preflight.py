"""One five-minute zero-solve Torch check before credit reoptimization."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

def main():
    assert sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID')
    import numpy as np
    start=time.monotonic()
    source=Path(__file__).parent
    out=Path(sys.argv[1]); out.mkdir(exist_ok=False)
    base=load(source/'run_fixed_price.py','credit_preflight_base')
    adapter=load(source/'natural_credit.py','credit_preflight_adapter')
    manifest,contract,objective,runtime,prepared,ref=base.authenticate(out)
    P,grid=ref['parameters'],np.asarray(ref['b_grid'])
    model=prepared.rt['model']
    q=adapter.construct(model,P,grid,float(ref['solution'].p_eq[0]))
    old_grid=grid.copy()
    if '--refine' in sys.argv:
        # Add economically derived continuation boundaries. Margins decrease
        # with age, guaranteeing resources for the next boundary plus the
        # unchanged kernel's 1e-6 positive-spending margin. Extra owner dust
        # puts liquidation strictly above the matching renter node.
        nodes=[float(q['minimum'][j].max()-q['sale'][h]+1e-4*(P.J-j)+1e-7*h)
               for j in range(P.J) for h in range(len(q['sale']))]
        grid=np.unique(np.r_[grid,[x for x in nodes if grid[0]<x<grid[-1]]])
        old_indices=np.searchsorted(grid,old_grid)
        np.testing.assert_array_equal(grid[old_indices],old_grid)
        q=adapter.construct(model,P,grid,float(ref['solution'].p_eq[0]))
        (out/'credit_grid.json').write_text(json.dumps(dict(grid=grid.tolist(),old_indices=old_indices.tolist(),
            method='Add worst-income natural beginning-wealth boundaries with explicit positive-spending margins',
            original_grid_preserved_exactly=True))+'\n')
    else:
        old_indices=np.arange(len(grid))
    # Independent general minimum-over-tenure feasibility recurrence.
    M=np.zeros((P.J,len(P.z_grid),len(q['cost'])))
    for j in range(P.J-1,-1,-1):
        L=np.empty(len(q['cost']))
        for h in range(len(L)):
            candidates=[-q['sale'][h]] if q['survival'][j]<1 else []
            if q['survival'][j]>0: candidates.append(M[j+1,:,h].max())
            L[h]=max(candidates)
        np.testing.assert_allclose(L,q['human'][j]-q['sale'],atol=2e-13,rtol=0)
        for zz in range(len(P.z_grid)):
            for old in range(len(L)):
                M[j,zz,old]=min((L[new]+q['oc'][new]-q['income'][j,zz])/P.R_gross
                    -(0 if new==old else q['sale'][old]-q['cost'][new]) for new in range(len(L)))
        np.testing.assert_allclose(M[j],q['minimum'][j,:,None]-q['sale'][None,:],atol=2e-13,rtol=0)
    # Twelve deterministic strict-support fixtures, including exact nodes.
    fixtures=[(-.1,False),(0.,False),(.25,False),(.5,True),(.75,True),(1.,True),(1.1,False),
              (.5,True),(.6,True),(.4,False),(.001,False),(.999,True)]
    for x,want in fixtures:
        got=adapter.reachable(np.array([0.,.5,1.]),np.array([False,True,True]),np.array([x]))[0]
        assert bool(got)==want
    shape=q['pre'].transpose(3,2,0,1)[:,:,None,:,:,None,None]
    old_g=np.asarray(ref['stationary_g_pre'])
    g=np.zeros((len(grid),)+old_g.shape[1:]); g[old_indices]=old_g
    bad=np.broadcast_to(~shape,g.shape)
    badmass=float(g[bad].sum())
    old_entry=np.asarray(P.fixed_reference_entry_conditional)
    entry=np.zeros((len(grid),old_entry.shape[1])); entry[old_indices]=old_entry
    conditional_bad=[float(entry[:,z][~q['pre'][0,z,0]].sum()) for z in range(len(P.z_grid))]
    economics_ok=grid[:,None,None,None] > (q['minimum'].T[None,None,:,:]-q['sale'][None,:,None,None])
    # economics_ok is wealth, tenure, income, age.
    economic_bad=np.broadcast_to(~economics_ok.transpose(0,1,3,2)[:,:,None,:,:,None,None],g.shape)
    economic_badmass=float(g[economic_bad].sum())
    ages=P.age_start+np.arange(P.J)*P.period_years
    result=dict(reference_label=base.LABEL,status='ready' if badmass<=2e-10 and max(conditional_bad)==0 else 'blocked_before_solve',
        slurm_job=os.environ['SLURM_JOB_ID'],model_solves=0,source_sha256={p.name:base.sha(p) for p in source.glob('*.py')},
        reference_manifest_sha256=base.MANIFEST_SHA,independent_general_recursion_pass=True,
        boolean_support_fixtures=12,baseline_mass_outside_discrete_solvency=badmass,
        baseline_mass_outside_economic_solvency=economic_badmass,
        infeasible_entry_conditional_mass_by_income=conditional_bad,
        rows=[dict(age=float(ages[j]),economic_renter_saving_floor=float(q['human'][j]),
             numerical_renter_saving_floor=float(q['floors'][j,0]),
             maximum_grid_tightening=float(np.max(q['floors'][j]-(q['human'][j]-q['sale'])))) for j in range(P.J)],
        grid=dict(nodes=len(grid),minimum=float(grid[0]),maximum=float(grid[-1])),
        elapsed_seconds=time.monotonic()-start)
    (out/'preflight.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__': main()
