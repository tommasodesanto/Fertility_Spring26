"""Recover a saved-path cross-section using finite household lifetimes.

No equilibrium solve. A 2023 household lives at most J more dated decisions,
so the exact saved prices through its last possible age suffice. Every recovered
2007--2023 aggregate row must reproduce the existing full-horizon replay.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'


def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def profile(e,P,grid):
    import numpy as np
    if P.child_state_mode!='independent_count': raise ValueError('Dependent-count observer required')
    g=e.g_current; rows=[]; large=[]
    weights=np.arange(P.n_parity,dtype=float);weights[-1]=P.tfr_top_bin_weight
    for j in range(P.J):
        x=g[:,:,:,j,:,:,:]; mass=float(x.sum())
        rooms=float(np.sum(x[:,0]*e.policy.hR_pol[:,0,:,j,:,:,:]))
        cap=float(np.sum(x[:,0]*np.minimum(e.policy.hR_pol[:,0,:,j,:,:,:],9)))
        big_all=big_children=0.
        for t,h in enumerate(P.H_own,1):
            amount=float(x[:,t].sum());rooms+=amount*h;cap+=amount*min(h,9)
            if h>=6:
                big_all+=amount
                big_children+=float(x[:,t,:,: ,1:,1:].sum())
        # Childless readiness states may use cs=1; they are not dependent children.
        children=float(x[:,:,:,:,1:,1:].sum())
        ceb=float(np.sum(x*weights[None,None,None,None,:,None]))/mass
        rows.append(dict(age=float(P.age_start+j*P.da),age_width=float(P.da),households=mass,
            owners=float(x[:,1:].sum()),rooms=rooms,capped_rooms=cap,with_children=children,
            children_ever_born=ceb,consumption=float(np.sum(x*e.policy.c_pol[:,:,:,j,:,:,:])),
            wealth=float(np.sum(x*grid[:,None,None,None,None,None]))))
        large.append(dict(age=float(P.age_start+j*P.da),age_width=float(P.da),
            without_children=big_all-big_children,with_children=big_children))
    totals={k:sum(r[k] for r in rows) for k in ('households','owners','rooms','capped_rooms','with_children','consumption','wealth')}
    if abs(totals['rooms']-float(e.demand_by_loc.sum()))>2e-10:
        raise ValueError('Age-profile housing does not aggregate to native demand')
    return dict(rows=rows,large_owner_age_cells=large,number_children_mass=g.sum(axis=(0,1,2,3,4,6)),
        totals=totals,child_observer='Positive dependent count among households with children ever born; excludes childless readiness states')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True);a=ap.parse_args()
    m=read(a.manifest)
    for p,h in m['file_sha256'].items():
        if sha(p)!=h:raise ValueError('Changed pinned input: '+p)
    import numpy as np
    sys.path.insert(0,m['replay_source'])
    import replay_e5f_original_queue_diagnostics as replay
    replay.np=np
    spec=read(m['spec']);sys.path.insert(0,str(Path(spec['batch'])/'source'))
    import run_e5f_original_queue_experiments as runner
    c=runner.load_context(m['spec'])
    old_smoke=read(c.spec['smoke_summary'])
    if old_smoke['status']!='passed' or old_smoke['spec_sha256']!=sha(m['spec']):
        raise ValueError('Frozen original-loop smoke must match')
    frozen,latest,prices,pensions,transfers=replay.validate_snapshot(Path(m['snapshot']))
    prices=np.asarray(prices,dtype=float)
    endpoint,reference=replay.endpoint_and_reference(Path(m['endpoint']),Path(m['endpoint_receipt']),Path(m['reference']),spec)
    import run_e5f_transition_calibration as fertility
    from e5f_balanced_terminal import _household_checks
    out=Path(m['output'])
    if out.exists():raise ValueError('Refusing to overwrite recovery')
    out.mkdir(parents=True)
    deadline=min(float(m['deadline_unix']),time.time()+1200);stop=threading.Event()
    count=5;backward_count=count+int(c.old.parameters.J)
    if backward_count!=22:raise ValueError('Expected 17 age cells and 2023 at date4')
    c.driver.save(out/'startup.json',dict(forward_dates=count,backward_dates=backward_count,
        finite_lifetime_argument='At date5, youngest age0 reaches last ageJ-1 at date5+J-1; last-age Bellman ignores continuation V. Thus saved V5 and all earlier dates are exact for the fixed100-date price path.',
        original_price_horizon=100,root_solves=0,deadline_unix=deadline,manifest_sha256=sha(a.manifest)))
    def heartbeat():
        while not stop.wait(60):
            c.driver.save(out/'heartbeat.json',dict(remaining_seconds=deadline-time.time()))
            if time.time()>=deadline:
                c.driver.save(out/'failure.json',dict(error='Recovery deadline'));os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    observations=[]; saved={}; pf=c.joined.pf
    def observe(i,e,P,grid,shared):
        rents=np.asarray(pf.rents_from_asset_prices(prices[:count],float(prices[count]),c.old.parameters),dtype=float)
        _,gates=_household_checks(e,P,shared,grid,float(rents[i]),c.primitive,c.audit)
        if not gates or not all(gates.values()):raise RuntimeError('Native household audit failed')
        observations.append(dict(period=i,calendar_year=2007+4*i,**fertility.period_fertility_diagnostics(e,P)))
        c.driver.save(out/'latest_date.json',dict(year=2007+4*i))
        if i==4:
            saved.update(calendar_year=2023,profile=profile(e,P,grid),fertility=observations[-1])
            c.driver.save(out/'model_2023.json',saved)
            with gzip.open(out/'native_2023_snapshot.pkl.gz','wb',compresslevel=1) as f:
                pickle.dump(dict(parameters=P,b_grid=grid,evaluation=e,shared=shared,
                    supply_rule=c.old.supply_rule,calendar_year=2023,source='permanent_shock_iteration3'),f,protocol=pickle.HIGHEST_PROTOCOL)
    try:
        with c.queue.original_queue_adapter(),c.cache.policy_cache(pf,max_bytes=12*1024**3):
            rents=np.asarray(pf.rents_from_asset_prices(prices[:backward_count],float(prices[backward_count]),c.old.parameters),dtype=float)
            calls=[0];native=pf.solve_date_policy
            def counted(**kwargs):
                result=native(**kwargs);calls[0]+=1
                c.driver.save(out/'backward_progress.json',dict(completed=calls[0],total=backward_count))
                return result
            with patch.object(pf,'solve_date_policy',counted):
                values,_=pf.backward_value_path(prices=prices[:backward_count],rents=rents,
                    psi_path=np.full(backward_count,spec['permanent_psi']),terminal_V=endpoint.policy.V,
                    base_parameters=c.old.parameters,b_grid=c.old.b_grid,transfer_path=transfers[:backward_count],
                    pension_path=pensions[:backward_count],payroll_tax_path=np.full(backward_count,.179))
            # Persist the continuation required for this exact recovery.
            with gzip.open(out/'continuation_2027.pkl.gz','wb',compresslevel=1) as f:
                pickle.dump(dict(V=values[count],price=float(prices[count])),f,protocol=pickle.HIGHEST_PROTOCOL)
            result=c.queue.queue_path(inherited=c.rebated.InheritedState(2007,c.old.initial_state),old_state=c.old,
                prices=prices[:count],pensions=pensions[:count],transfers=transfers[:count],psi=spec['permanent_psi'],
                terminal=NS(parameters=NS(psi_child=spec['permanent_psi']),policy=NS(V=values[count]),asset_price=float(prices[count])),observer=observe)
        check=replay.compare_rows(result.rows,frozen[:count]);check['rows']=count
        if len(observations)!=count or saved.get('calendar_year')!=2023:raise RuntimeError('Missing requested cross-section')
        c.driver.save(out/'rows.json',result.rows);c.driver.save(out/'fertility.json',observations)
        c.driver.save(out/'verification.json',dict(status='PASS',observed_year=2023,row_reproduction=check,
            backward_dates=backward_count,forward_dates=count,model_snapshot_sha256=sha(out/'native_2023_snapshot.pkl.gz'),
            model_source_sha256=sha(out/'model_2023.json'),finite_equilibrium_converged=False,root_solves=0,
            source='Same fixed-price permanent-shock iteration3 as the mock transition graphs'))
        c.driver.verify_pins(m['file_sha256'])
    except BaseException as exc:
        c.driver.save(out/'failure.json',dict(error_type=type(exc).__name__,error=str(exc)));raise
    finally:stop.set()


if __name__=='__main__':main()
