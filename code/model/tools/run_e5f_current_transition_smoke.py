#!/usr/bin/env python3
"""Bounded native two-date transition mapping, not an equilibrium path solve."""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'
os.environ['MPLBACKEND'] = 'Agg'
import argparse
import contextlib
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import signal
import shutil
import subprocess
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TOOLS = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def queue_values(queue):
    return np.asarray(queue.due_in_16 + queue.due_in_20 if hasattr(queue, 'due_in_16') else queue)


def exact_difference(first, second):
    result = {}
    for name in ('prices','rents'):
        result[name] = float(np.max(np.abs(getattr(first,name)-getattr(second,name))))
    if len(first.values) != len(second.values) or len(first.rows) != len(second.rows):
        raise ValueError('Different cache replay dimensions')
    result['values'] = max(float(np.max(np.abs(x-y))) for x,y in zip(first.values,second.values))
    result['terminal_g'] = float(np.max(np.abs(first.terminal_state.g_pre-second.terminal_state.g_pre)))
    for name in ('scheduled_entries','scheduled_raw_entries'):
        result[name] = float(np.max(np.abs(queue_values(getattr(first.terminal_state,name))-queue_values(getattr(second.terminal_state,name)))))
    gaps = []
    for left,right in zip(first.rows,second.rows):
        if left.keys() != right.keys(): raise ValueError('Different replay row fields')
        for key,value in left.items():
            if isinstance(value,(int,float,np.number)):
                if not np.isfinite(value) or not np.isfinite(right[key]):raise ValueError('Nonfinite dated result')
                gaps.append(abs(float(value)-float(right[key])))
            elif value != right[key]:raise ValueError('Different replay metadata: '+key)
    result['rows'] = max(gaps,default=0.)
    return result


def trial_prices(reference_price, terminal_price, arm):
    if not np.isfinite([reference_price,terminal_price]).all() or min(reference_price,terminal_price)<=0:
        raise ValueError('Positive finite endpoint prices required')
    if arm == 'reference':return np.full(2,reference_price)
    if arm != 'natural' or reference_price == terminal_price:
        raise ValueError('Natural diagnostic requires genuinely different endpoint price')
    # Interior geometric points avoid a single-date jump; no equilibrium claim.
    return np.exp(np.log(reference_price)+np.array([1/3,2/3])*np.log(terminal_price/reference_price))


def audit_gates(audit):
    budget,purchase,estate = audit['budget'],audit['purchase'],audit['estate']
    gates = {'budget':float(budget['budget_excess_mass']) <= 2e-10,
             'estate_funded':estate['status']=='funded',
             'estate_next_cohort':estate['audit_id']=='estate_funded_dated_entry_provisional_net_v1',
             'negative_estates':float(estate['estate']['totals']['net_negative']) <= 1e-10}
    for key in ('end_mortgage_floor_violation_mass','purchase_threshold_violation_mass','transaction_outside_grid_mass',
                'negative_estate_exposure_mass','saving_outside_grid_mass'):
        if key in purchase:gates[key]=float(purchase[key]) <= 2e-10
    gates['transaction_wealth']=float(purchase['maximum_occupied_transaction_wealth_error']) <= 1e-9
    return gates


def collect_pins(reference_case,terminal_case,contract,approval,terminal_approval):
    approved=read(approval); terminal=read(terminal_approval)
    paths={Path(p) for p in approved['current_source_files']}
    paths.update(Path(p) for p in terminal['current_source_files'])
    paths.update(TOOLS/name for name in ('run_e5f_current_transition_smoke.py',
        'e5f_current_transition_runtime.py','e5f_exact_policy_cache.py',
        'e5f_overnight_estate_audit.py','e5f_solvency_credit_benchmark.py'))
    paths.update([Path(contract),Path(approval),Path(terminal_approval),
        Path(approved['baseline_receipt_path']),Path(terminal['terminal_receipt_path'])])
    for case in (reference_case,terminal_case):
        paths.update([Path(case)/'initial_state.pkl.gz',Path(case)/'receipt.json',Path(case)/'parameters.csv'])
    return {str(p.resolve()):sha(p) for p in sorted(paths)}


def archive_sources(output,pins):
    manifest={}
    for source,digest in pins.items():
        p=Path(source)
        if p.suffix!='.py':continue
        relative=p.relative_to(ROOT) if p.is_relative_to(ROOT) else Path('external')/digest/p.name
        target=Path(output)/'source_snapshot'/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        if sha(target)!=digest:raise ValueError('Source changed during snapshot: '+source)
        target.chmod(0o444);manifest[source]=dict(snapshot=str(target.resolve()),sha256=digest)
    write(Path(output)/'source_snapshot_manifest.json',manifest)


def verify_pins(pins):
    for path,expected in pins.items():
        if sha(path)!=expected:raise ValueError('Pinned source/input changed: '+path)


def validate_approval(path, reference_checkpoint):
    approval=read(path)
    if approval.get('approved') is not True:raise ValueError('Lead-reviewed baseline approval required')
    if approval.get('reference_checkpoint_sha256')!=sha(reference_checkpoint):raise ValueError('Approved reference differs')
    receipt=Path(approval['baseline_receipt_path'])
    if sha(receipt)!=approval['baseline_receipt_sha256']:raise ValueError('Approved baseline receipt differs')
    status=read(receipt).get('status')
    if status not in ('native_reference_replay_pass','native_reference_tables_pass_array_review_required'):
        raise ValueError('Native baseline tables have not passed')
    if status=='native_reference_tables_pass_array_review_required' and approval.get('arrays_reviewed') is not True:
        raise ValueError('Explicit lead array review required')
    pins=approval.get('current_source_files',{})
    if not pins:raise ValueError('Missing approved current-source pins')
    verify_pins(pins)
    return approval


def validate_terminal_approval(path,terminal_case,baseline_approval):
    approval=read(path)
    if approval.get('approved') is not True or approval.get('status')!='native_terminal_replay_verified':
        raise ValueError('Explicit verified native terminal approval required')
    if approval.get('terminal_checkpoint_sha256')!=sha(Path(terminal_case)/'initial_state.pkl.gz'):
        raise ValueError('Approved native terminal checkpoint differs')
    receipt=Path(approval['terminal_receipt_path'])
    if sha(receipt)!=approval['terminal_receipt_sha256']:raise ValueError('Approved terminal receipt differs')
    if receipt.resolve()!=(Path(terminal_case)/'receipt.json').resolve():raise ValueError('Terminal receipt path mismatch')
    pins=approval.get('current_source_files',{})
    baseline=read(baseline_approval)['current_source_files']
    model_pins={p:v for p,v in baseline.items() if '/intergen_eqscale_seq_optimized/' in p}
    if not model_pins or any(pins.get(p)!=v for p,v in model_pins.items()):
        raise ValueError('Native terminal and baseline model source pins differ')
    verify_pins(pins)
    return approval


def prepare(args):
    for name in ('reference_case','terminal_case','contract','approval','terminal_approval'):
        if getattr(args,name) is None:raise ValueError('Required --'+name.replace('_','-'))
    validate_approval(args.approval,args.reference_case/'initial_state.pkl.gz')
    validate_terminal_approval(args.terminal_approval,args.terminal_case,args.approval)
    args.output.mkdir(parents=True,exist_ok=False)
    plan=dict(source_pins=collect_pins(args.reference_case,args.terminal_case,args.contract,args.approval,args.terminal_approval),
        reference_case=str(args.reference_case.resolve()),terminal_case=str(args.terminal_case.resolve()),
        contract=str(args.contract.resolve()),approval=str(args.approval.resolve()),terminal_approval=str(args.terminal_approval.resolve()),deadline=min(time.time()+1200,args.deadline or float('inf')),
        case_seconds=600,workers=1,arms=['reference','natural'],horizon_dates=2,cache_max_bytes=2*1024**3,
        expected_solve_count='16 dated Bellman requests; exact cache can reduce to12 actual calls',
        scope='Two-date prescribed-price operator/cache test, not solved transition; natural prices change across dates',
        unchanged=['psi_child','payroll_tax','housing_supply','initial_g_pre','grid','entry_distribution'],
        economic_change='Native solvency credit enabled only in natural arm; inherited gross donor utility and net estate funding')
    archive_sources(args.output,plan['source_pins'])
    write(args.output/'plan.json',plan)


def evaluate(args):
    plan=read(args.output.parent/'plan.json');verify_pins(plan['source_pins'])
    if time.time()>=plan['deadline']:raise TimeoutError('Smoke deadline expired')
    args.output.mkdir(parents=True,exist_ok=False)
    import e5f_current_transition_runtime as native
    import e5f_exact_policy_cache as cache
    validate_approval(plan['approval'],Path(plan['reference_case'])/'initial_state.pkl.gz')
    validate_terminal_approval(plan['terminal_approval'],Path(plan['terminal_case']),plan['approval'])
    loaded=native.setup(args.output/'runtime',contract=Path(plan['contract']),reference=Path(plan['reference_case']))
    packet,rt=loaded['selected'],loaded['runtime']
    pf=rt['primitive'].pf
    terminal=packet
    if args.arm=='natural':
        tables=[]
        for key in ('reference_case','terminal_case'):
            with (Path(plan[key])/'parameters.csv').open() as stream:
                tables.append({r['parameter']:float(r['estimate']) for r in csv.DictReader(stream)})
        if len(tables[0]) != 31 or tables[0] != tables[1]:raise ValueError('All 31 terminal/reference parameters must match exactly')
        case=Path(plan['terminal_case']);receipt=read(case/'receipt.json')
        if sha(case/'initial_state.pkl.gz')!=receipt['case_checkpoint_sha256']:raise ValueError('Terminal checkpoint differs')
        with gzip.open(case/'initial_state.pkl.gz','rb') as stream:terminal=pickle.load(stream)
    P=copy.deepcopy(loaded['parameters']);P.native_solvency_credit=args.arm=='natural'
    for name in ('psi_child','tau_pay','H0','r_bar','xi_supply','property_tax_lump_sum_transfer'):
        if not np.array_equal(np.asarray(getattr(P,name)),np.asarray(getattr(terminal['parameters'],name))):
            raise ValueError('Terminal changes retained primitive: '+name)
    if not np.array_equal(packet['b_grid'],terminal['b_grid']):raise ValueError('Terminal grid differs')
    g0=packet['stationary_g_pre'];grid=packet['b_grid'];E=float(g0[:,:,:,0,:,:,:].sum())
    initial=pf.stationary_initial_state(g0,E,float(packet['evaluation'].births),P,1/2.1)
    if not np.array_equal(initial.g_pre,g0):raise ValueError('Initial distribution changed')
    p0=float(packet['solution'].p_eq[0]);pT=float(terminal['solution'].p_eq[0])
    prices=trial_prices(p0,pT,args.arm);rents=pf.rents_from_asset_prices(prices,pT,P)
    kwargs=dict(prices=prices,psi_path=np.full(2,P.psi_child),terminal_price=pT,
        terminal_V=terminal['evaluation'].policy.V,base_parameters=P,b_grid=grid,initial_state=initial,
        supply_rule=packet['supply_rule'],birth_to_entry_conversion=1/2.1,
        transfer_path=np.full(2,P.property_tax_lump_sum_transfer),pension_path=np.full(2,P.pension),
        payroll_tax_path=np.full(2,P.tau_pay))
    results=[];timings=[];audits=[];stats=None
    import e5f_overnight_estate_audit as estate
    import e5f_solvency_credit_benchmark as credit_audit
    for enabled in (False,True):
        arm_audits=[]
        def observer(period,e,parameters,b_grid,shared,next_entrant_cohort):
            if time.time()>=plan['deadline']:raise TimeoutError('Dated mapping deadline')
            budget=rt['primitive'].dated_budget(e,parameters,shared,b_grid,float(rents[period]))
            purchase=(credit_audit.audit_purchase_accounting if args.arm=='natural' else rt['accounting'].audit_purchase_accounting)(e,parameters,shared,b_grid,rt['model'])
            try:funding=estate.audit(e,parameters,b_grid,next_entrant_cohort=next_entrant_cohort)
            except estate.EstateFundingShortfall as exc:funding=exc.audit
            row=dict(period=period,cache_enabled=enabled,budget=budget,purchase=purchase,estate=funding)
            row['gates']=audit_gates(row);arm_audits.append(row)
            write(args.output/('audits_cache_on.json' if enabled else 'audits_cache_off.json'),arm_audits)
            write(args.output/'progress.json',dict(epoch=time.time(),period=period,cache_enabled=enabled))
            if not all(row['gates'].values()):raise RuntimeError('Dated household/accounting gates fail')
        start=time.monotonic()
        with (cache.policy_cache(pf,max_bytes=plan['cache_max_bytes']) if enabled else contextlib.nullcontext()) as cache_stats:
            result=pf.evaluate_path_at_prices(**kwargs,dated_observer=observer)
            if enabled:stats=cache_stats.snapshot()
        timings.append(time.monotonic()-start);results.append(result);audits.append(arm_audits)
        write(args.output/('path_cache_on.json' if enabled else 'path_cache_off.json'),result.rows)
        write(args.output/'timing.json',dict(seconds=timings,cache=stats))
    clean_audits=[[{k:v for k,v in row.items() if k!='cache_enabled'} for row in group] for group in audits]
    if clean_audits[0]!=clean_audits[1]:raise RuntimeError('Cache replay dated audits differ')
    differences=exact_difference(*results)
    write(args.output/'replay_comparison.json',differences)
    if any(x!=0 for x in differences.values()):raise RuntimeError('Exact cache replay differs')
    result=results[-1]
    gates=dict(mass=result.maximum_mass_accounting_error<2e-8,
               policy_reproduction=result.maximum_policy_reproduction_error<1e-10,
               projection=result.maximum_feasibility_projection_mass<=1e-6)
    fiscal=max(abs(float(r['scaled_pension_budget_residual'])) for r in result.rows)
    if args.arm=='reference':gates.update(market=result.maximum_market_residual<=2e-4,fiscal=fiscal<=1e-6)
    write(args.output/'gates.json',gates)
    if not all(gates.values()):raise RuntimeError('Native mapping gates fail')
    with (args.output/'path.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=result.rows[0]);writer.writeheader();writer.writerows(result.rows)
    verify_pins(plan['source_pins'])
    write(args.output/'complete.json',dict(status='passed_native_mapping_cache_smoke',arm=args.arm,
        scope=plan['scope'],prices=prices.tolist(),terminal_price=pT,cache=stats,seconds=timings,
        exact_differences=differences,gates=gates,maximum_market_residual=result.maximum_market_residual,
        maximum_fiscal_residual=fiscal,equilibrium_transition=False,initial_distribution_sha256=hashlib.sha256(g0.tobytes()).hexdigest(),
        fixed_psi_child=P.psi_child,fixed_payroll_tax=P.tau_pay,source_pins=plan['source_pins']))


def run(args):
    plan=read(args.output/'plan.json');verify_pins(plan['source_pins']);records=[]
    if (args.output/'checkpoint.json').exists():raise ValueError('Do not overwrite prior run')
    for arm in plan['arms']:
        if time.time()>=plan['deadline']:break
        with (args.output/(arm+'.log')).open('x') as log:
            proc=subprocess.Popen([sys.executable,__file__,'--stage','evaluate','--arm',arm,'--output',str(args.output/arm)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            deadline=min(plan['deadline'],time.time()+plan['case_seconds'])
            while proc.poll() is None and time.time()<deadline:
                write(args.output/'heartbeat.json',dict(epoch=time.time(),arm=arm,pid=proc.pid,completed=len(records)));time.sleep(5)
            timed_out=proc.poll() is None
            if timed_out:
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait(timeout=5)
        ok=proc.returncode==0 and (args.output/arm/'complete.json').exists()
        records.append(dict(arm=arm,returncode=proc.returncode,status='passed' if ok else 'timeout' if timed_out else 'failed'))
        write(args.output/'checkpoint.json',dict(records=records))
        if not ok:break
    write(args.output/'complete.json',dict(status='passed' if len(records)==2 and all(r['status']=='passed' for r in records) else 'failed',records=records))


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['prepare','run','evaluate'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--arm',choices=['reference','natural'])
    p.add_argument('--terminal-approval',type=Path);p.add_argument('--approval',type=Path);p.add_argument('--contract',type=Path);p.add_argument('--reference-case',type=Path);p.add_argument('--terminal-case',type=Path);p.add_argument('--deadline',type=float)
    args=p.parse_args();{'prepare':prepare,'run':run,'evaluate':evaluate}[args.stage](args)

if __name__=='__main__':main()
