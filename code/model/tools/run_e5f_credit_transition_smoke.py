#!/usr/bin/env python3
"""Bounded fixed-parameter transition mapping and exact-cache regression.

This is a two-date operator/speed smoke, not a solved equilibrium transition.
Economic changes are restricted to the approved credit/solvency experiment.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
os.environ['MPLBACKEND']='Agg'
import argparse, contextlib, copy, csv, gzip, hashlib, json, pickle, signal, subprocess, sys, time
from pathlib import Path
import numpy as np
from run_e5f_credit_benchmark import module, read, sha, write
ROOT=Path(__file__).resolve().parents[3]
REFERENCE=ROOT/'tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/search/de_0093/case'
CONTRACT=REFERENCE.parents[2]/'production_contract.json'
BENCHMARK=ROOT/'output/model/daytime_calibration_20260927/credit_benchmark/run_v2/natural/case'


def exact_difference(a,b):
    differences={}
    for name in ('prices','rents'):
        differences[name]=float(np.max(np.abs(getattr(a,name)-getattr(b,name))))
    differences['values']=max(float(np.max(np.abs(x-y))) for x,y in zip(a.values,b.values))
    differences['terminal_g']=float(np.max(np.abs(a.terminal_state.g_pre-b.terminal_state.g_pre)))
    for name in ('scheduled_entries','scheduled_raw_entries'):
        differences[name]=float(np.max(np.abs(np.asarray(getattr(a.terminal_state,name))-getattr(b.terminal_state,name))))
    assert len(a.rows)==len(b.rows)
    gaps=[]
    for x,y in zip(a.rows,b.rows):
        assert x.keys()==y.keys()
        for k,v in x.items():
            if isinstance(v,(int,float)):
                gaps.append(abs(float(v)-float(y[k])))
            else: assert v==y[k],k
    differences['rows']=max(gaps,default=0.)
    return differences


def evaluate(args):
    plan=read(args.output.parent/'plan.json')
    for path,expected in plan['source_pins'].items():assert sha(path)==expected,path
    c=read(CONTRACT)
    os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=sha(CONTRACT)
    os.environ['E5F_LOCAL_EXECUTION_AUTHORIZATION']=c['execution']['authorization_id']
    from inspect_e5f_saved_households import load
    packet,rt,reference_receipt=load(CONTRACT,REFERENCE,args.output)
    pf=rt['primitive'].pf
    assert pf.calendar.model is rt['model']
    import e5f_credit_transition_runtime as adapter
    import e5f_solvency_credit_benchmark as credit
    import e5f_exact_policy_cache as cache
    queue_install=adapter.install_split_queue(pf,args.output/'runtime')
    initial=adapter.initial_state(pf,packet)
    assert np.array_equal(initial.g_pre,packet['stationary_g_pre'])
    terminal=packet
    if args.arm=='natural':
        assert sha(BENCHMARK/'initial_state.pkl.gz')==read(BENCHMARK/'receipt.json')['case_checkpoint_sha256']
        with gzip.open(BENCHMARK/'initial_state.pkl.gz','rb') as stream:terminal=pickle.load(stream)
        credit.install(rt['model'],args.output/'runtime',enabled=True)
        adapter.install_continuation(rt['model'],args.output/'runtime')
    P=copy.deepcopy(packet['parameters']);grid=packet['b_grid'];H=2
    pterminal=float(terminal['solution'].p_eq[0]);price=np.full(H,pterminal)
    # This constant-price boundary is only for operator/cache verification.
    kwargs=dict(prices=price,psi_path=np.full(H,P.psi_child),terminal_price=pterminal,
        terminal_V=terminal['evaluation'].policy.V,base_parameters=P,b_grid=grid,
        initial_state=initial,supply_rule=packet['supply_rule'],birth_to_entry_conversion=1/2.1,
        transfer_path=np.full(H,P.property_tax_lump_sum_transfer),
        pension_path=np.full(H,P.pension),payroll_tax_path=np.full(H,P.tau_pay))
    original=pf.calendar.evaluate_period;audits=[]
    estate=module('transition_estate_audit',c['files']['estate_audit']['path'])
    def observed(*a,**kw):
        e=original(*a,**kw);parameters=a[2];shared=a[4]
        rent=float(pf.rents_from_asset_prices(price,pterminal,parameters)[len(audits)%H])
        budget=rt['primitive'].dated_budget(e,parameters,shared,grid,rent)
        purchase=(credit.audit_purchase_accounting if args.arm=='natural' else rt['accounting'].audit_purchase_accounting)(e,parameters,shared,grid,rt['model'])
        try:funding=estate.audit(e,parameters,grid)
        except estate.EstateFundingShortfall as exc:funding=exc.audit
        audits.append(dict(budget=budget,purchase=purchase,estate=funding))
        write(args.output/'progress.json',dict(epoch=time.time(),arm=args.arm,forward_dates=len(audits),total_forward_dates=4))
        return e
    pf.calendar.evaluate_period=observed
    results=[];times=[];stats=None
    try:
        for enabled in (False,True):
            start=time.monotonic()
            with (cache.policy_cache(pf,max_bytes=2*1024**3) if enabled else contextlib.nullcontext()) as current:
                result=pf.evaluate_path_at_prices(**kwargs)
                if enabled:stats=current.snapshot()
            results.append(result);times.append(time.monotonic()-start)
    finally:pf.calendar.evaluate_period=original
    write(args.output/'raw_mapping_rows.json',[r.rows for r in results])
    write(args.output/'raw_timing.json',dict(seconds=times,cache=stats))
    write(args.output/'audits.json',audits)
    differences=exact_difference(*results)
    assert all(v==0 for v in differences.values()),differences
    for result in results:
        assert result.maximum_mass_accounting_error<2e-8
        assert result.maximum_policy_reproduction_error<1e-10
        assert result.maximum_feasibility_projection_mass<=1e-6
    if args.arm=='reference':
        assert results[0].maximum_market_residual<=2e-4
        assert max(abs(r['scaled_pension_budget_residual']) for r in results[0].rows)<=1e-6
    rows=results[-1].rows
    with (args.output/'path.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    write(args.output/'audits.json',audits)
    result=dict(status='passed_operator_and_exact_cache_smoke',arm=args.arm,horizon_dates=H,
        scope='Two-date fixed-price mapping only; not equilibrium transition or terminal certification',
        cache_off_seconds=times[0],cache_on_seconds=times[1],speed_ratio=times[0]/times[1],cache=stats,
        exact_differences=differences,maximum_market_residual=results[-1].maximum_market_residual,
        queue=queue_install,parameter_changes=[],fixed_psi_child=P.psi_child,
        initial_distribution_sha256=hashlib.sha256(initial.g_pre.tobytes()).hexdigest(),
        terminal_boundary='reference' if args.arm=='reference' else 'conditional fixed-entry benchmark for operator smoke only',
        estate_price_timing='inherited current decision-price after saving; no timing reform',
        source_pins=plan['source_pins'])
    write(args.output/'complete.json',result)


def run(args):
    args.output.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__),Path(__file__).with_name('e5f_credit_transition_runtime.py'),Path(__file__).with_name('e5f_solvency_credit_benchmark.py'),Path(__file__).with_name('e5f_exact_policy_cache.py'),Path(__file__).with_name('run_e5f_credit_benchmark.py'),Path(__file__).with_name('inspect_e5f_saved_households.py')]
    plan=dict(source_pins={str(p.resolve()):sha(p) for p in sources},deadline=min(time.time()+1200,args.deadline or float('inf')),
        case_seconds=600,workers=1,arms=['reference','natural'],horizon_dates=2,
        expected_solve_count='16 dated Bellman requests; exact cache may reduce to12 actual calls',
        economic_scope='Approved credit experiment only; current split16/20 queue, same pre-shock distribution, no historical reweighting')
    write(args.output/'plan.json',plan);records=[]
    for arm in plan['arms']:
        with (args.output/(arm+'.log')).open('w') as log:
            proc=subprocess.Popen([sys.executable,__file__,'--stage','evaluate','--arm',arm,'--output',str(args.output/arm)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            deadline=min(plan['deadline'],time.time()+plan['case_seconds'])
            while proc.poll() is None and time.time()<deadline:
                write(args.output/'heartbeat.json',dict(epoch=time.time(),arm=arm,pid=proc.pid,completed=len(records)));time.sleep(5)
            if proc.poll() is None:
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait(timeout=5)
        ok=proc.returncode==0 and (args.output/arm/'complete.json').exists()
        records.append(dict(arm=arm,returncode=proc.returncode,status='passed' if ok else 'failed'))
        write(args.output/'checkpoint.json',dict(records=records))
        if not ok:break
    write(args.output/'complete.json',dict(status='passed' if len(records)==2 and all(r['status']=='passed' for r in records) else 'failed',records=records))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['run','evaluate'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--arm',choices=['reference','natural']);p.add_argument('--deadline',type=float);args=p.parse_args()
    run(args) if args.stage=='run' else evaluate(args)
