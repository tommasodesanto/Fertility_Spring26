"""Bounded staged recalibration using the previously verified case adapter.

Two exact-reference histories smoke-test this exact execution/collection loop
before any candidate search. No model, objective, or candidate rule lives here.
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

import run_e5f_bounded_calibration_refinement as adapter
import build_e5f_bounded_refinement_plan as planner


def run_case_process(command, log, timeout, active, lock):
    with open(log, 'w') as stream:
        p = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        with lock:
            active.add(p)
        try:
            try:
                code = p.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                return {'returncode': 124, 'timeout': True}
            return {'returncode': code, 'timeout': False}
        finally:
            with lock:
                active.discard(p)


def execute(contract_path, expected):
    adapter.verify(contract_path, expected)
    c = adapter.read_json(contract_path)
    for path, sha in c['code_sha256'].items():
        adapter.verify(path, sha)
    if c['max_workers'] != 23 or c['max_cases'] != 39 or c['case_timeout_seconds'] != 6000 or c['total_seconds'] != 28800:
        raise RuntimeError('Undeclared search budget')
    base_path = Path(c['base_plan'])
    base = adapter.load_plan(base_path, c['base_plan_sha256'])
    if base.get('choice_model') != 'fertility_nest':
        raise RuntimeError('Wrong experimental model')
    for path, sha in c['reference_sha256'].items():
        adapter.verify(path, sha)
    root = Path(c['output'])
    if base_path.parent != root/'smoke' or {p.name for p in root.iterdir()} != {'smoke'}:
        raise RuntimeError('Search output must contain only the predeclared smoke inputs')
    start = time.time()
    deadline = min(start+c['total_seconds'], float(base['launch_deadline_epoch']))
    active, lock, completed = set(), threading.Lock(), []
    def stop(signum, frame):
        with lock:
            for p in list(active):
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        raise RuntimeError(f'Controller interrupted by signal {signum}')
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    count = 0
    def stage(plan_path):
        nonlocal count
        plan_sha = adapter.digest(plan_path)
        plan = adapter.load_plan(plan_path, plan_sha)
        if time.time()+c['case_timeout_seconds']+300 > deadline:
            raise RuntimeError('Insufficient reserved time for a complete stage')
        count += len(plan['cases'])
        if count > c['max_cases']:
            raise RuntimeError('Case budget exceeded')
        failure = []
        with cf.ThreadPoolExecutor(max_workers=c['max_workers']) as pool:
            futures = {}
            for case in plan['cases']:
                command = [sys.executable, str(Path(adapter.__file__).resolve()), '--plan', str(plan_path),
                           '--plan-sha256', plan_sha, '--case-id', str(case['id'])]
                future = pool.submit(run_case_process, command, plan_path.parent/f"case_{case['id']:03d}.log",
                                     c['case_timeout_seconds'], active, lock)
                futures[future] = case
            pending = set(futures)
            while pending:
                done, pending = cf.wait(pending, timeout=60, return_when=cf.FIRST_COMPLETED)
                for future in done:
                    case = futures[future]
                    result = future.result()
                    row = dict(stage=plan['stage'], case_id=case['id'], label=case['label'], **result)
                    out = plan_path.parent/case['output']
                    if result['returncode'] == 0:
                        receipt = adapter.read_json(out/'case_receipt.json')
                        if receipt['status'] != 'complete' or receipt['plan_sha256'] != plan_sha:
                            raise RuntimeError('Invalid completed-case receipt')
                        for name in ('summary.json','target_fit_long.csv','parameter_table.csv'):
                            adapter.verify(out/name, receipt['artifact_sha256'][name])
                        adapter.validate_result(out, plan, case)
                        row.update(loss=receipt['loss'], summary=str(out/'summary.json'))
                        adapter.write_json(root/'latest_completed_case.json',row)
                        best_path = root/'live_best.json'
                        best = adapter.read_json(best_path) if best_path.exists() else {'loss': float('inf')}
                        if row['loss'] < best['loss']:
                            adapter.write_json(best_path,row)
                    else:
                        failure.append(row)
                    completed.append(row)
                    adapter.write_json(root/'completed_cases.json',completed)
                adapter.write_json(root/'heartbeat.json',dict(stage=plan['stage'],elapsed_seconds=time.time()-start,
                    epoch=time.time(),completed_cases=len(completed),running_cases=len(pending),failed_cases=len(failure)))
                for future in pending:
                    case = futures[future]
                    h = plan_path.parent/case['output']/'heartbeat.json'
                    if h.exists() and time.time()-h.stat().st_mtime > 1800:
                        stop(signal.SIGTERM,None)
        if failure:
            raise RuntimeError(f'Stage stopped after failed cases: {failure}')
        return planner.collect(plan_path, plan_sha, require_complete=True)
    try:
        # The reference smoke uses the same subprocess, receipts and collector
        # as the later search stages. Its exact-reference gates stay mandatory.
        prior,status,rows = stage(base_path)
        for action, folder in ((planner.coordinate,'coordinates'),(planner.joint,'joint'),(planner.repeats,'repeats')):
            dest = root/folder
            prior = dict(prior,parent_plan_sha256=status['plan_sha256'])
            action(prior,status,rows,dest)
            prior,status,rows = stage(dest/'plan.json')
        if len(rows)!=2 or rows[0]['loss']!=rows[1]['loss']:
            raise RuntimeError('Final exact repeats disagree')
        selected = adapter.read_json(root/'best_so_far.json')
        selected_dir = Path(selected['best']['summary']).parent
        reference_dir = Path(c['comparison_reference'])
        reference_summary = adapter.read_json(reference_dir/'summary.json')
        selected_summary = adapter.read_json(selected_dir/'summary.json')
        if reference_summary['target_fingerprint'] != selected_summary['target_fingerprint']:
            raise RuntimeError('Comparison changes target contract')
        combined = []
        selected_rows = adapter.read_csv(selected_dir/'target_fit_long.csv')
        reference_rows = adapter.read_csv(reference_dir/'target_fit_long.csv')
        if len(selected_rows) != 12 or len(reference_rows) != 12:
            raise RuntimeError('Incomplete comparison target table')
        for new, old in zip(selected_rows,reference_rows):
            if any(new[k] != old[k] for k in ('moment','target','weight')):
                raise RuntimeError('Comparison changes targets or weights')
            combined.append(dict(moment=new['moment'], target=new['target'], weight=new['weight'],
                sequential_model=old['model'], nested_model=new['model'],
                sequential_gap=old['gap'], nested_gap=new['gap'],
                sequential_loss=old['loss_contribution'],nested_loss=new['loss_contribution']))
        planner.write_csv(root/'comparison_target_fits.csv',combined)
        parameter_rows=[]
        for name,folder in (('sequential',reference_dir),('selected_nested',selected_dir)):
            parameter_rows.extend(dict(model=name,**row) for row in adapter.read_csv(folder/'parameter_table.csv'))
        planner.write_csv(root/'comparison_parameters.csv',parameter_rows)
        adapter.write_json(root/'final_receipt.json',dict(status='complete',elapsed_seconds=time.time()-start,
            completed_cases=len(completed),selected=selected,
            sequential_reference_loss=reference_summary['best_candidate']['transition_loss'],
            relative_loss_change=selected['loss']/reference_summary['best_candidate']['transition_loss']-1,
            repetitions=rows,contract_sha256=expected,production_promoted=False))
    except BaseException as error:
        adapter.write_json(root/'failure.json',dict(error=repr(error),elapsed_seconds=time.time()-start,
                                                  completed_cases=len(completed)))
        raise

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contract',type=Path,required=True)
    p.add_argument('--contract-sha256',required=True)
    a=p.parse_args()
    execute(a.contract.resolve(),a.contract_sha256)
