"""Torch-only two-arm price-start diagnosis; no normalization or model changes.

The unchanged PAYGO wrapper certifies/rejects the actual solver return. A proxy
saves that return's diagnostics first, so a strict-gate failure retains evidence.
This is not a calibrated point, production solution, or model-existence test.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback
from types import SimpleNamespace


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def plain(x):
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    if hasattr(x,'tolist'):return plain(x.tolist())
    if isinstance(x,float) and not math.isfinite(x):return {'nonfinite':repr(x)}
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    return repr(x)


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(plain(value),indent=2,allow_nan=False)+'\n');temp.replace(path)


def authenticate(args):
    if sys.platform!='linux' or not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Torch Slurm execution only; no local numerical execution')
    if digest(args.plan)!=args.plan_sha256:raise RuntimeError('Plan pin differs')
    plan=json.loads(args.plan.read_text())
    if plan['status']!='lead_approved':raise RuntimeError('Plan awaits explicit lead review')
    if plan['case_id']!='initial_0041_block' or plan['arms']!=[1.,.5]:raise RuntimeError('Wrong matched case/arms')
    if not 0<plan['case_cap_seconds']<=1800 or not 0<plan['parent_cap_seconds']<=3600:
        raise RuntimeError('Budget enlarged')
    for key in ('runner','contract','checkpoint'):
        pin=plan[key]
        if not pin['sha256'] or digest(pin['path'])!=pin['sha256']:raise RuntimeError('Bad pin: '+key)
    if Path(plan['runner']['path']).resolve()!=Path(__file__).resolve():raise RuntimeError('Wrong wrapper')
    if time.time()>=float(plan['global_end_epoch']):raise RuntimeError('Author deadline expired')
    contract=json.loads(Path(plan['contract']['path']).read_text())
    records=json.loads(Path(plan['checkpoint']['path']).read_text())['records']
    matches=[r for r in records if r['case']==plan['case_id']]
    if len(matches)!=1:raise RuntimeError('Case identity ambiguous')
    record=matches[0]
    if record['status']!='inadmissible' or record['design']!='broad_coverage' or record['lane']!='block':
        raise RuntimeError('Original failed-case classification changed')
    if record['error']['error']!='Initial housing equilibrium failed its unchanged strict gate':
        raise RuntimeError('Original failure is not the matched gate')
    if 'initial = evaluate' not in record['error']['traceback']:raise RuntimeError('Not first-intercept failure')
    if record['error']['context']['contract_sha256']!=plan['contract']['sha256']:
        raise RuntimeError('Diagnostic contract is not the original failed-case contract')
    return plan,contract,record


def child(args):
    plan,contract,record=authenticate(args)
    if args.arm not in (0,1):raise RuntimeError('Invalid arm')
    if time.time()>=args.deadline:raise RuntimeError('Case deadline expired before setup')
    if args.deadline>min(time.time()+plan['case_cap_seconds'],float(plan['global_end_epoch'])):
        raise RuntimeError('Child deadline exceeds approved cap')
    output=args.output;output.mkdir(parents=True,exist_ok=False)
    write(output/'progress.json',dict(stage='authenticating_native_runtime',epoch=time.time()))
    # Heavy imports occur only in this supervised Torch child.
    import numpy as np
    import e5f_evening_calibration_runtime as runtime
    contract=dict(contract,objective=contract['lanes']['block']['objective'])
    objective=json.loads(Path(contract['objective']['path']).read_text())
    evaluator=runtime.setup(contract,objective,output)
    evaluator.case_output=output;evaluator.due=True
    P=evaluator.bind(record['point'])
    P.psi_child=float(contract['normalization']['initial_psi'])
    if P.psi_child!=0.14281100340255604:raise RuntimeError('Original fixed intercept changed')
    original=np.asarray(evaluator.selected['solution'].p_eq,dtype=float)
    factor=plan['arms'][args.arm];start_price=original*factor
    if len(original)!=1 or not np.isfinite(start_price).all():raise RuntimeError('Bad price start')
    tax,fiscal_rule=evaluator.adapter.pension_tax_from_demographics(P)
    write(output/'request.json',dict(case=record['case'],point=record['point'],fixed_psi=P.psi_child,
        original_start=original,price_start=start_price,factor=factor,due=True,
        calibrated=False,normalization_performed=False,deadline=args.deadline,plan_sha256=args.plan_sha256))
    write(output/'target_definition.json',objective)
    # Preserve named original-reference parameters, explicitly not this failed point's fit.
    reference=Path(contract['reference_case'])/'parameters.csv'
    if reference.exists():(output/'original_reference_parameters.csv').write_bytes(reference.read_bytes())
    write(output/'proposed_parameters.json',dict(point=record['point'],fixed_psi=P.psi_child,
        public_parameters={k:v for k,v in vars(P).items() if not k.startswith('_')}))
    real=evaluator.rt['model'].solve_markov_income_equilibrium
    calls=0
    def capture(*pos,**kw):
        nonlocal calls
        calls+=1
        if calls!=1:raise RuntimeError('One GE call per arm only')
        write(output/'progress.json',dict(stage='solving_fixed_intercept_GE',epoch=time.time()))
        begun=time.monotonic()
        sol,params,price=real(*pos,**kw)
        receipt=dict(status='solver_return_not_certification',seconds=time.monotonic()-begun,
            converged=bool(sol.converged),price=price,tolerance=params.tol_eq,
            timings=sol.timings,mean_children=float(sol.mean_parity),
            completed_fertility=2*float(sol.mean_parity),own_rate=float(sol.own_rate),
            calibrated=False,full_scientific_audit_performed=False)
        write(output/'solver_return.json',receipt)
        write(output/'returned_parameters.json',{k:v for k,v in vars(params).items() if not k.startswith('_')})
        write(output/'progress.json',dict(stage='unchanged_PAYGO_strict_gate',epoch=time.time()))
        return sol,params,price
    try:
        sol,params,price,fiscal=evaluator.rt['solve_balanced_initial_equilibrium'](
            model=SimpleNamespace(solve_markov_income_equilibrium=capture),parameters=P,
            b_grid=evaluator.selected['b_grid'],initial_prices=start_price,payroll_tax=tax,
            marginal_tolerance=1e-9,fiscal_tolerance=1e-6,warm_price_state={})
        evaluator.adapter.verify_pension_ratio(fiscal)
        # Identical native gate plus explicit finite evidence; never certify NaN.
        if not math.isfinite(float(sol.timings['best_eq_error'])):raise RuntimeError('Nonfinite returned residual')
        runtime.verify_sources(contract)
        write(output/'result.json',dict(status='housing_and_PAYGO_gates_passed_only',fiscal=fiscal,
            ge_calls=calls,calibrated=False,full_scientific_audit_performed=False))
    except Exception as exc:
        write(output/'result.json',dict(status='rejected_or_error',error_type=type(exc).__name__,error=str(exc),
            traceback=traceback.format_exc(),ge_calls=calls,calibrated=False))
        raise


def parent(args):
    plan,contract,record=authenticate(args)
    args.output.mkdir(parents=True,exist_ok=False)
    start=time.time();end=min(start+plan['parent_cap_seconds'],float(plan['global_end_epoch']))
    receipts=[]
    for arm in (0,1):
        now=time.time();deadline=min(now+plan['case_cap_seconds'],end-10.)
        if deadline<=now:
            receipts.append(dict(arm=arm,status='unrun_parent_budget'));continue
        out=args.output/f'arm_{arm}'
        cmd=[sys.executable,str(Path(__file__).resolve()),'--plan',str(args.plan.resolve()),
             '--plan-sha256',args.plan_sha256,'--output',str(out.resolve()),
             '--arm',str(arm),'--deadline',str(deadline)]
        env=dict(os.environ)
        for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS'):env[k]='1'
        with (args.output/f'arm_{arm}.log').open('w') as log:
            proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            timed_out=False
            while proc.poll() is None:
                remaining=deadline-time.time()
                if remaining<=0:
                    timed_out=True
                    try:os.killpg(proc.pid,signal.SIGTERM)
                    except ProcessLookupError:pass
                    try:proc.wait(timeout=2.)
                    except subprocess.TimeoutExpired:pass
                    # Kill the owned group even if its leader exited, removing descendants.
                    try:os.killpg(proc.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    proc.wait();break
                write(args.output/'heartbeat.json',dict(arm=arm,pid=proc.pid,epoch=time.time(),deadline=deadline))
                try:proc.wait(timeout=min(30.,remaining))
                except subprocess.TimeoutExpired:pass
        receipts.append(dict(arm=arm,factor=plan['arms'][arm],status='censored_timeout' if timed_out else 'child_finished',
            returncode=proc.returncode,elapsed=time.time()-now,deadline=deadline))
        write(args.output/'checkpoint.json',receipts)
    write(args.output/'complete.json',dict(scope='two-arm fixed-intercept numerical diagnosis',
        cases=receipts,elapsed=time.time()-start,calibrated=False,plan_sha256=args.plan_sha256))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha256',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--arm',type=int)
    p.add_argument('--deadline',type=float,default=0.)
    a=p.parse_args()
    child(a) if a.arm is not None else parent(a)
if __name__=='__main__':main()
