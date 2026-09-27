#!/usr/bin/env python3
"""Fixed-parameter credit experiment on the authenticated overnight solution.

The child-benefit level is held fixed. Birth renewal and estate entry funding
are reported, not repaired with parameter changes. Household numerical gates
remain active. Every subprocess has its own bounded lifetime and output folder.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import time

import numpy as np


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temp.replace(path)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def evaluate(args):
    root = args.output.parent
    plan = read(root / 'plan.json')
    assert sha(__file__) == plan['runner_sha256']
    assert sha(plan['helper']) == plan['helper_sha256']
    contract = Path(plan['contract'])
    assert sha(contract) == plan['contract_sha256']
    c = read(contract)
    os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256'] = sha(contract)
    os.environ['E5F_LOCAL_EXECUTION_AUTHORIZATION'] = c['execution']['authorization_id']
    driver = module('credit_benchmark_frozen_driver', c['files']['driver']['path'])
    c, objective = driver.verify(contract)
    driver.verify_execution(c)
    reference_case = Path(plan['reference_case'])
    receipt = read(reference_case / 'receipt.json')
    assert sha(reference_case / 'initial_state.pkl.gz') == receipt['case_checkpoint_sha256']
    args.output.mkdir(exist_ok=False)
    (args.output / 'runtime').mkdir()
    runtime, tax, selected, rt, *_ = driver.setup(c, objective, receipt['point'], args.output / 'runtime')
    with gzip.open(reference_case / 'initial_state.pkl.gz', 'rb') as stream:
        reference = pickle.load(stream)
    helper = module('credit_solvency_adapter', plan['helper'])
    installation = helper.install(rt['model'], args.output / 'runtime', enabled=args.arm == 'natural')
    write(args.output / 'installation.json', installation)

    # Keep estate shortfalls visible while retaining complete benchmark tables.
    original_estate_audit = runtime.estate.audit
    def estate_with_report(*a, **kw):
        try:
            return original_estate_audit(*a, **kw)
        except runtime.estate.EstateFundingShortfall as exc:
            return exc.audit
    runtime.estate.audit = estate_with_report
    if args.arm == 'natural':
        rt['accounting'].audit_purchase_accounting = helper.audit_purchase_accounting

    def fixed_solve(point, output, deadline):
        P = copy.deepcopy(reference['parameters'])
        fixed_psi = float(P.psi_child)
        payroll_tax, fiscal_rule = runtime.adapter.pension_tax_from_demographics(P)
        start = time.monotonic()
        write(output / 'stationary_solves.json', [dict(status='started', epoch=time.time(), psi_child=fixed_psi)])
        sol, P, price, fiscal = rt['solve_balanced_initial_equilibrium'](
            model=rt['model'], parameters=P, b_grid=reference['b_grid'],
            initial_prices=reference['solution'].p_eq, payroll_tax=payroll_tax,
            marginal_tolerance=1e-9, fiscal_tolerance=1e-6, warm_price_state={})
        assert float(P.psi_child) == fixed_psi
        runtime.adapter.verify_pension_ratio(fiscal)
        seconds = time.monotonic() - start
        fertility = float(rt['chain'].extract_moments(sol, P)['tfr'])
        normalization = dict(status='fixed_child_benefit_no_renormalization', psi_child=fixed_psi,
            target=2.1, completed_fertility=fertility, absolute_gap=abs(fertility-2.1),
            stationary_solves=1, stationary_solve_seconds=seconds)
        ledger = dict(status='completed', seconds=seconds, psi_child=fixed_psi,
            price=float(price[0]), market_residual=float(sol.timings['best_eq_error']),
            price_evaluations=int(sol.timings['unique_fast_price_evaluations']))
        write(output / 'stationary_solves.json', [ledger])
        renewal = dict(status='reported_not_imposed_in_fixed_parameter_comparison',
            entry_rate=float(sol.entry_rate), adjusted_birth_children=float(sol.adult_entry_adjusted_birth_children),
            relative_replacement_gap=float(sol.adult_entry_stationary_relative_gap),
            closure='Conditional stationary cross-section with inherited normalized entry; no closed-renewal claim')
        return (sol, P, price, seconds, normalization), fiscal_rule, renewal

    runtime.normalize = fixed_solve
    result = runtime.evaluate_point(tax=tax, objective=objective, selected=selected, runtime=rt,
        point=receipt['point'], output=args.output / 'case', deadline_epoch=plan['deadline'], graphs=True)
    result.update(status='fixed_parameter_credit_diagnostic', arm=args.arm,
        benchmark_plan_sha256=sha(root / 'plan.json'), adapter_sha256=plan['helper_sha256'],
        normalization_inputs={'mode':'fixed_reference_psi_child','psi_child':float(reference['parameters'].psi_child)},
        economic_changes=[] if args.arm=='reference' else [
            'Replace purchase/collateral/unsecured limits by grid-resolved continuation feasibility and net-estate solvency at every possible death date.'],
        scientific_scope='Same-parameter stationary credit diagnostic; not a recalibration or a fully closed demographic counterfactual',
        approximation='Conservative feasible-grid boundary; native value cutoff may reject extremely negative finite utility')
    with gzip.open(args.output / 'case/initial_state.pkl.gz','rb') as stream:
        packet = pickle.load(stream)
    totals = result['estate_funding']['estate']['totals']
    if args.arm == 'natural':
        assert totals['net_negative'] <= 1e-9, totals
    P = packet['parameters']
    reference_values = vars(reference['parameters'])
    changed = []
    for key, value in vars(P).items():
        if key.startswith('_') or key not in reference_values:
            continue
        old = reference_values[key]
        try:
            equal = np.array_equal(value, old, equal_nan=True)
        except (TypeError, ValueError):
            equal = repr(value) == repr(old)
        if not equal:
            changed.append(key)
    result['parameter_fields_changed_from_reference'] = changed
    # eq_iter is overwritten by the equilibrium solver on each price evaluation.
    economic_changed = [key for key in changed if key != 'eq_iter']
    assert not economic_changed, 'Same-parameter benchmark changed fields: ' + repr(economic_changed)
    write(args.output / 'case/receipt.json', result)
    # Free-coordinate estimates remain fixed; update descriptive labels only.
    path = args.output / 'case/parameters.csv'
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if row['parameter'] in receipt['point'] or row['parameter']=='psi_child':
            row['status']='fixed at reference; no recalibration or fertility normalization'
        if row['parameter']=='child_benefit_CRRA_coefficient':
            row['status']='derived from fixed reference one-child benefit'
        if args.arm=='natural' and row['parameter']=='financed_share':
            row['status']='retained parameter value; collateral rule inactive in solvency benchmark'
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    write(args.output / 'complete.json', result)


def run(args):
    out=args.output;out.mkdir(parents=True,exist_ok=False)
    original=read(args.contract);c=copy.deepcopy(original)
    ref=read(args.reference_case/'receipt.json')
    c['initial_point']=ref['point'];c['execution']['authorization_id']='tommaso_daytime_solvency_credit_20260927'
    c['budget']['absolute_end_epoch']=min(time.time()+1800, args.deadline or float('inf'))
    c['files']['credit_benchmark_runner']={'path':str(Path(__file__).resolve()),'sha256':sha(__file__)}
    c['files']['credit_benchmark_helper']={'path':str(args.helper.resolve()),'sha256':sha(args.helper)}
    c['files']['credit_benchmark_reference_receipt']={'path':str((args.reference_case/'receipt.json').resolve()),'sha256':sha(args.reference_case/'receipt.json')}
    contract=out/'contract.json';write(contract,c)
    plan=dict(contract=str(contract.resolve()),contract_sha256=sha(contract),runner_sha256=sha(__file__),
        helper=str(args.helper.resolve()),helper_sha256=sha(args.helper),reference_case=str(args.reference_case.resolve()),
        deadline=c['budget']['absolute_end_epoch'],case_seconds=600,workers=1,
        cases=['reference','natural','natural_repeat'],fixed_child_benefit=True,
        author_authority='September27: proceed with same-parameter credit benchmark and explicit solvency at death')
    write(out/'plan.json',plan)
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[name]='1'
    env['MPLBACKEND']='Agg'
    records=[]
    for name in plan['cases']:
        if time.time()>=plan['deadline']:break
        arm='reference' if name=='reference' else 'natural'
        with (out/(name+'.log')).open('w') as log:
            p=subprocess.Popen([sys.executable,__file__,'--stage','evaluate','--output',str(out/name),'--arm',arm],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            limit=min(plan['deadline'],time.time()+plan['case_seconds'])
            while p.poll() is None and time.time()<limit:
                write(out/'heartbeat.json',dict(epoch=time.time(),case=name,pid=p.pid,completed=len(records),deadline=limit))
                time.sleep(5)
            timed_out=p.poll() is None
            if timed_out:
                os.killpg(p.pid,signal.SIGTERM)
                try:p.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait(timeout=5)
        success=p.returncode==0 and not timed_out and (out/name/'complete.json').exists()
        record=dict(case=name,status='success' if success else ('timeout' if timed_out else 'failed'),returncode=p.returncode)
        if success:
            result=read(out/name/'complete.json');record['loss']=result['loss']
            if name=='reference':
                with (args.reference_case/'target_fit.csv').open() as f:old={r['moment']:r for r in csv.DictReader(f)}
                with (out/name/'case/target_fit.csv').open() as f:new={r['moment']:r for r in csv.DictReader(f)}
                gaps={k:float(new[k]['model'])-float(v['model']) for k,v in old.items()}
                write(out/'reference_comparison.json',dict(model_differences=gaps,tolerance=1e-5))
                with (args.reference_case/'parameters.csv').open() as f: old_parameters={r['parameter']:float(r['estimate']) for r in csv.DictReader(f)}
                with (out/name/'case/parameters.csv').open() as f: new_parameters={r['parameter']:float(r['estimate']) for r in csv.DictReader(f)}
                assert old_parameters == new_parameters, 'Parameter replay differs'
                success=max(abs(v) for v in gaps.values())<=1e-5
                record['status']='success' if success else 'reference_reproduction_failed'
            elif name=='natural_repeat':
                for filename in ('target_fit.csv','parameters.csv'):
                    if (out/'natural/case'/filename).read_bytes()!=(out/name/'case'/filename).read_bytes():
                        success=False;record['status']='repeat_mismatch'
        records.append(record);write(out/'checkpoint.json',dict(records=records))
        if not success:break
    write(out/'complete.json',dict(status='complete' if len(records)==3 and all(r['status']=='success' for r in records) else 'incomplete',records=records,epoch=time.time()))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['run','evaluate'],required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--contract',type=Path)
    parser.add_argument('--reference-case',type=Path);parser.add_argument('--helper',type=Path)
    parser.add_argument('--arm',choices=['reference','natural'])
    parser.add_argument('--deadline',type=float,help='Preserve an earlier absolute experiment deadline on a numerical repair')
    args=parser.parse_args()
    if args.stage=='run':
        if not all((args.contract,args.reference_case,args.helper)):parser.error('run requires contract, reference-case and helper')
        run(args)
    else:evaluate(args)
