"""Explicit horizon diagnostic conditional on a fitted six-period 2023 state.

No historical refit, structural change, or production promotion. The native
two-policy six-period replay must pass before separate longer forecasts launch.
"""
import argparse
import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import threading
import time

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np

def extend_coordinates(values, old_count, count, shift=0):
    a=np.asarray(values,dtype=float)
    if a.shape!=(3*(old_count+1),) or not np.isfinite(a).all() or np.any(a<=0):
        raise ValueError('Expected complete positive price/pension/rebate blocks')
    if not 0<=shift<=old_count or count<1:raise ValueError('Invalid numerical seed horizon')
    indices=np.minimum(np.arange(count+1)+shift,old_count)
    return a.reshape(3,old_count+1)[:,indices].reshape(-1).copy()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--spec',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--mode',choices=['smoke','forecast'],required=True)
    ap.add_argument('--policy',choices=['baseline_rebate','tax2_rebate'],default='baseline_rebate')
    ap.add_argument('--seed',choices=['flat','long'],default='flat');args=ap.parse_args()
    read=lambda p:json.loads(Path(p).read_text())
    spec=read(args.spec);manifest=read(spec['runtime_manifest']);plan=read(manifest['prior_plan'])
    sys.path[:0]=[spec['helper'],str(Path(plan['source_root'])/'code/model/tools'),str(Path(plan['source_root'])/'code/model'),spec['cache_directory']]
    import run_e5f_final_rebated_history as driver
    import e5f_rebated_surprises as rebated
    import e5f_exact_policy_cache as cache
    from e5f_rebated_initial_bridge import build_rebated_initial_state
    from e5f_balanced_terminal import TerminalAuditControls
    from run_e5f_successive_surprises_overnight import standard_graphs
    assert Path(driver.__file__).parent.resolve()==Path(spec['helper']).resolve()
    driver.verify_pins(spec['file_sha256']);driver.verify_pins(manifest['file_sha256']);driver.verify_pins(plan['file_sha256'])
    proof=read(spec['boundary_verification']);assert proof['status']=='passed'
    cache_proof=read(spec['cache_proof']);digest=driver.sha(cache.__file__)
    assert digest==cache_proof['cache_sha256']==spec['cache_sha256']
    assert cache_proof['status']=='verified' and cache_proof['exact_mapping_equal'] and cache_proof['household_and_accounting_mapping_valid']
    source=Path(spec['source_case']);contract=read(source/'contract_receipt.json')
    assert contract['count']==6 and contract['case']=='A0' and contract['manifest_sha256']==driver.sha(spec['source_manifest'])
    driver.verify_pins(read(spec['source_manifest'])['file_sha256'])
    fits=read(source/'realized_fit.json');assert [r['year'] for r in fits]==[2007,2011,2015,2019]
    assert read(source/'finite_history_complete.json')['realized']==fits and all(abs(r['gap'])<=.005 for r in fits)
    _,joined,*_=rebated._runtime()
    summary=read(manifest['initial_summary']);assert summary['status']=='verified_rebated_initial_smoke'
    item=summary['checkpoint'];cp=item.get('path',item.get('checkpoint'));expected=item.get('sha256',item.get('checkpoint_sha256'))
    assert driver.sha(cp)==expected==contract['initial_checkpoint_sha256']
    for pair in read(manifest['kernel_equivalence'])['pairs']:
        assert driver.sha(pair['initial'])==driver.sha(pair['history'])==pair['sha256']
    with gzip.open(cp,'rb') as f:packet=pickle.load(f)
    with gzip.open(source/'realized_state_2023.pkl.gz','rb') as f:inherited=pickle.load(f)
    assert inherited.year==2023
    initial_hash=hashlib.sha256(inherited.households.g_pre.tobytes()).hexdigest()
    raw=read(spec['initial_raw_summary'])
    old=build_rebated_initial_state(packet=packet,normalization=raw['normalization'],outside_origin_entry_share=plan['outside_origin_entry_share'],preference_change_2023=0.)
    demographics=driver.migration_case(packet['demographic_seed'],'A0')
    controls=dict(plan['history_root_controls']);controls.update(manifest.get('root_controls',{}))
    controls.setdefault('transfer_bounds',[1e-10,10.]);controls['max_evaluations']=4 if args.mode=='smoke' else min(24,int(controls['max_evaluations']))
    audit=TerminalAuditControls(**plan['terminal_template']['audit_controls'])
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    if (out/'contract_receipt.json').exists():raise ValueError('Refusing to overwrite an existing diagnostic')
    count=6 if args.mode=='smoke' else 24
    if args.mode=='forecast':
        smoke=read(spec['smoke_receipt']);assert smoke['status']=='passed' and smoke['spec_sha256']==driver.sha(args.spec)
    remaining=min(1200 if args.mode=='smoke' else 7200,spec['numerical_deadline_unix']-time.time())
    if remaining<120:raise TimeoutError('Insufficient time before the author cutoff')
    deadline=time.monotonic()+remaining;stop=threading.Event()
    receipt=dict(source_history=str(source),source_history_count=6,forecast_count=count,source_history_refitted=False,
        initial_g_pre_sha256=initial_hash,psi=fits[-1]['psi'],case='A0',start_year=2023,seed=args.seed,mode=args.mode,
        spec_sha256=driver.sha(args.spec),horizon_verified=False,production_eligible=False,results={})
    driver.save(out/'contract_receipt.json',receipt)
    def heartbeat():
        while not stop.wait(60):driver.save(out/'heartbeat.json',dict(remaining_seconds=max(0,deadline-time.monotonic())))
    threading.Thread(target=heartbeat,daemon=True).start()
    try:
        with cache.policy_cache(joined.pf,max_bytes=6*1024**3) as stats:
            for policy in (['baseline_rebate','tax2_rebate'] if args.mode=='smoke' else [args.policy]):
                tax=.01 if policy=='baseline_rebate' else .02
                original=source/'policies'/policy;root=read(original/'root_receipt.json')
                assert root['converged'] and root['finite_horizon_market_fiscal_converged'] and root['final_reproduction_max_abs']<=2e-10
                assert root['count']==6 and root['start_year']==2023 and root['psi']==fits[-1]['psi']
                initial=extend_coordinates(root['final']['prices'],6,count)
                if args.seed=='long' and count==24:
                    long_root=read(spec['long_seed_root']);assert long_root['count']==24 and long_root['start_year']==2007 and long_root['converged'] and long_root['finite_horizon_market_fiscal_converged']
                    base=read(source/'policies/baseline_rebate/root_receipt.json')
                    ratio=initial/extend_coordinates(base['final']['prices'],6,count)
                    initial=extend_coordinates(long_root['final']['prices'],24,count,shift=4)*ratio
                policy_old=copy.copy(old);policy_old.parameters=copy.deepcopy(old.parameters)
                P=policy_old.parameters;P.tau_H=4*tax;P.user_cost_rate=float(P.R_gross)+float(P.delta)+float(P.tau_H)-1.
                folder=out/policy
                driver.save(folder/'initial_coordinates.json',dict(prices=initial,seed=args.seed,seed_is_numerical_only=True))
                result,detail=driver.solve_forecast(inherited=inherited,old=policy_old,demographics=demographics,
                    psi=fits[-1]['psi'],count=count,initial=initial,controls=controls,audit=audit,deadline=deadline,folder=folder,case='A0')
                if result.next_state is None:raise RuntimeError('Conditional policy forecast did not converge')
                np.testing.assert_array_equal(detail['snapshot']['evaluation'].g_pre,inherited.households.g_pre)
                assert hashlib.sha256(inherited.households.g_pre.tobytes()).hexdigest()==initial_hash
                if args.mode=='smoke':
                    np.testing.assert_array_equal(result.root_receipt['final']['prices'],root['final']['prices'])
                    np.testing.assert_array_equal(result.root_receipt['final']['residual'],root['final']['residual'])
                    assert driver.clean(result.path.rows)==read(original/'rows.json')
                    assert driver.clean(detail['observations'])==read(original/'fertility.json')
                standard_graphs(detail['snapshot'],result,folder/'graphs')
                receipt['results'][policy]=dict(finite_converged=True,annual_tax=tax,same_initial_g_pre=True,root_sha256=driver.sha(folder/'root_receipt.json'),exact_original_replay=args.mode=='smoke')
                driver.save(folder/'summary.json',receipt['results'][policy]);driver.save(out/'latest_completed.json',receipt)
            receipt.update(status='passed',cache=stats.snapshot());driver.verify_pins(spec['file_sha256'])
            driver.save(out/'summary.json',receipt)
    except BaseException as exc:
        driver.save(out/'failure.json',dict(error_type=type(exc).__name__,error=str(exc)));raise
    finally:stop.set()

if __name__=='__main__':main()
