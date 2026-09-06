#!/usr/bin/env python3
"""Correct the reporting guard for an externally fixed first-child term.

The original preparation and numerical plans remain immutable. The fixed term
is checked in the fingerprinted summary, where the scientific driver stores it.
"""
import argparse,math,json
from pathlib import Path
import empirical_room_gap as run
adapter=run.adapter
original_new_plan=run.new_plan

def new_plan(*args):
    p=original_new_plan(*args)
    p['empirical_gap_collector_sha256']=adapter.digest(__file__)
    return p
run.new_plan=new_plan

def collect(path,sha,complete=True):
    plan=adapter.load_plan(path,sha)
    adapter.verify(run.__file__,plan['empirical_gap_preparation_sha256'])
    if 'empirical_gap_collector_sha256' in plan:
        adapter.verify(__file__,plan['empirical_gap_collector_sha256'])
    if len(plan['held_nine'])!=9:raise RuntimeError('Nine held coordinates required')
    for case in plan['cases']:
        out=path.parent/case['output'];receipt=out/'case_receipt.json'
        if not receipt.exists():continue
        r=adapter.read_json(receipt)
        for name in ('parameter_table.csv','summary.json'):
            adapter.verify(out/name,r['artifact_sha256'][name])
        values={v['parameter']:float(v['value']) for v in adapter.read_csv(out/'parameter_table.csv')}
        for key,value in plan['held_nine'].items():
            if not math.isclose(values[key],value,rel_tol=0,abs_tol=2e-12):raise RuntimeError(f'Held coordinate changed: {key}')
        adapter.verify(path.parent/case['center'],case['center_sha256'])
        center=adapter.read_json(path.parent/case['center'])['best_candidate']['theta']
        theta=adapter.read_json(out/'summary.json')['best_candidate']['theta']
        for key in ('hbar_child_rooms','hbar_first_child_jump'):
            if not math.isclose(theta[key],center[key],rel_tol=0,abs_tol=2e-12):raise RuntimeError(f'Proposed housing coordinate changed: {key}')
    result=run.profile.collect(path,sha,complete)
    adapter.write_json(path.parent/'report/collector_version.json',dict(
        collector_sha256=adapter.digest(__file__),preparation_sha256=plan['empirical_gap_preparation_sha256'],
        correction='Read the externally fixed first-child term from summary.best_candidate.theta, not the free-parameter table; all physical-coordinate checks and numerical gates unchanged'))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('collect','repeat'))
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--sha256',required=True);p.add_argument('--allow-partial',action='store_true');a=p.parse_args()
    plan,status,rows=collect(a.plan,a.sha256,not a.allow_partial)
    if a.action=='repeat':run.repeat(plan,rows)
    else:print(json.dumps(status))
