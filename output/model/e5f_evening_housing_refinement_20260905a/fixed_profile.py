#!/usr/bin/env python3
"""Prepare/collect a conditional bound diagnostic only after the bounded search.

No model is solved by preparation or collection. Numerical jobs use the
unchanged scientific driver with its existing explicit fixed-jump option.
"""
import argparse
import copy
import datetime
import json
import math
import shutil
import sys
import time
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path[:0] = [str(OUT),str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import fixed_jump_adapter as fixed
import run_e5f_bounded_calibration_refinement as base
import build_e5f_bounded_refinement_plan as planner


def collect(path, sha, complete=True):
    plan=fixed.load_plan(path,sha)
    fixed.verify(__file__,plan['profile_preparation_sha256'])
    rows,fits,parameters,missing,contracts=[],[],[],[],set()
    for case in plan['cases']:
        directory=path.parent/case['output']
        receipt_path=directory/'case_receipt.json'
        if not receipt_path.exists():missing.append(case['id']);continue
        receipt=fixed.read_json(receipt_path)
        if receipt['status']!='complete' or receipt['plan_sha256']!=sha or receipt['case_id']!=case['id']:
            raise RuntimeError('Profile receipt mismatch')
        for name,expected in receipt['artifact_sha256'].items():fixed.verify(directory/name,expected)
        s,models,gates=fixed.validate_result(directory,plan,case)
        if 'held_coordinates' in plan:
            measured={r['parameter']:float(r['value']) for r in fixed.read_csv(directory/'parameter_table.csv')}
            if len(plan['held_coordinates'])!=9:raise RuntimeError('Expected nine unchanged coordinates')
            for name,value in plan['held_coordinates'].items():
                if not math.isclose(measured[name],value,rel_tol=0,abs_tol=2e-12):
                    raise RuntimeError(f'Conditional profile changed held coordinate {name}')
            theta=s['best_candidate']['theta']
            if abs(theta['hbar_first_child_jump']+3*theta['hbar_child_rooms']-plan['preserved_three_child_floor'])>2e-12:
                raise RuntimeError('Conditional profile changed the three-child floor')
        import collect_e5f_transition_calibration as validator
        shared={k:copy.deepcopy(s[k]) for k in ('target_measurements','dated_measurement_contract','renewal_accounting_contract','external_closure_contract','model_profile')}
        # The only varying external profile setting is individually checked above.
        del shared['model_profile']['first_child_room_jump']
        shared['population_bridge']=validator.population_bridge_contract(s['population_bridge'])
        contracts.add(json.dumps(shared,sort_keys=True))
        rows.append(dict(id=case['id'],label=case['label'],loss=s['best_candidate']['transition_loss'],
            jump=case['fixed_first_child_room_jump'],first_birth_rooms=models['housing_increment_0to1'],
            mean_rooms=models['aggregate_mean_occupied_rooms_18_85'],
            summary=str((directory/'summary.json').resolve()),summary_sha256=fixed.digest(directory/'summary.json'),
            plan=str(path.resolve()),plan_sha256=sha,elapsed_seconds=receipt['elapsed_seconds']))
        fits.extend(dict(case_id=case['id'],label=case['label'],**r) for r in fixed.read_csv(directory/'target_fit_long.csv'))
        parameters.extend(dict(case_id=case['id'],label=case['label'],**r) for r in fixed.read_csv(directory/'parameter_table.csv'))
    if len(contracts)>1:raise RuntimeError('Mixed profile contracts beyond the declared jump')
    report=path.parent/'report';report.mkdir(exist_ok=True)
    rows.sort(key=lambda r:(r['loss'],r['id']))
    for name,data in (('all_candidates',rows),('all_target_fits',fits),('all_parameters',parameters)):
        planner.write_csv(report/f'{name}.csv',data)
    status=dict(status='partial' if missing else 'complete',completed_cases=len(rows),missing_cases=missing,
        plan_sha256=sha,best=rows[0] if rows else None,scope='conditional_profile_not_calibration',production_promoted=False)
    fixed.write_json(report/'summary.json',status)
    if rows:
        fixed.write_json(OUT/'profile_latest_completed.json',dict(stage=plan['stage'],latest=max(rows,key=lambda r:Path(r['summary']).stat().st_mtime)))
        best=OUT/'profile_best_so_far.json'
        if not best.exists() or rows[0]['loss']<fixed.read_json(best)['best']['loss']:
            fixed.write_json(best,dict(stage=plan['stage'],best=rows[0],scope='conditional_profile_not_calibration',production_promoted=False))
    if complete and missing:raise RuntimeError(f'Incomplete profile stage: {missing}')
    return plan,status,rows


def new_plan(prior,out,stage):
    plan=planner.fresh_plan(prior,out,stage)
    plan.update(schema='e5f_fixed_first_child_profile_v1',scope='conditional_profile_not_calibration',
        adapter_sha256=fixed.digest(fixed.__file__),profile_preparation_sha256=fixed.digest(__file__),
        launch_deadline_epoch=datetime.datetime(2026,9,6,0,20,tzinfo=datetime.timezone.utc).timestamp())
    return plan


def copy_reference(source,out):
    s=fixed.read_json(source);ref=out/'reference';ref.mkdir()
    for name in ('summary.json','target_fit_long.csv','parameter_table.csv'):
        shutil.copy2(source.parent/name,ref/name)
    name=s['best_candidate']['candidate'];(ref/'cases'/name).mkdir(parents=True)
    shutil.copy2(source.parent/'cases'/name/'transition_path.csv',ref/'cases'/name/'transition_path.csv')
    return ref/'summary.json'


def repeated_plan(prior,best,out,stage):
    chosen_path=Path(best['plan']);chosen=fixed.read_json(chosen_path)
    base.verify(chosen_path,best['plan_sha256'])
    case=next(c for c in chosen['cases'] if c['id']==best['id'])
    center_path=chosen_path.parent/case['center'];base.verify(center_path,case['center_sha256'])
    reference_path=Path(best['summary']);base.verify(reference_path,best['summary_sha256'])
    jump=fixed.read_json(reference_path)['best_candidate']['theta']['hbar_first_child_jump']
    plan=new_plan(prior,out,stage);reference=copy_reference(reference_path,out)
    plan['input_sha256']={str(p.relative_to(out)):fixed.digest(p) for p in reference.parent.rglob('*') if p.is_file()}
    for i in range(2):
        planner.add_case(plan,out,fixed.read_json(center_path),f'exact_fixed_jump_repeat_{i+1}',radius=case['radius'])
        plan['cases'][-1].update(fixed_first_child_room_jump=jump,reference=str(reference.relative_to(out)),reference_sha256=fixed.digest(reference))
    planner.finish_plan(plan,out)


def smoke():
    latest_start=datetime.datetime(2026,9,5,23,50,tzinfo=datetime.timezone.utc).timestamp()
    if time.time()>latest_start:raise RuntimeError('Conditional-profile start window closed')
    path=OUT/'repeats/plan.json'
    if not path.exists():raise RuntimeError('Complete unchanged-bound search first')
    prior,status,rows=planner.collect(path,base.digest(path))
    if len(rows)!=2:raise RuntimeError('Require both bounded final repeats')
    for row in rows:
        receipt=base.read_json(Path(row['summary']).parent/'case_receipt.json')
        if not receipt['reference']['exact_twelve_row_fit']:raise RuntimeError('Unverified bounded winner')
    best=base.read_json(OUT/'best_so_far.json')['best']
    s=base.read_json(best['summary'])
    if not math.isclose(s['best_candidate']['theta']['hbar_first_child_jump'],.5,rel_tol=0,abs_tol=1e-12):
        raise RuntimeError('Contingent ceiling hypothesis not met')
    repeated_plan(prior,best,OUT/'profile_smoke','conditional_profile_two_exact_smokes')


def grid(prior,rows):
    if prior['stage']!='conditional_profile_two_exact_smokes' or len(rows)!=2:
        raise RuntimeError('Two profile smokes required')
    for row in rows:
        if not fixed.read_json(Path(row['summary']).parent/'case_receipt.json')['reference']['exact_twelve_row_fit']:
            raise RuntimeError('Missing exact profile smoke')
    out=OUT/'profile_grid';plan=new_plan(prior,out,'conditional_profile_three_child_floor_preserved')
    s=fixed.read_json(rows[0]['summary']);theta=s['best_candidate']['theta']
    j0,c0=theta['hbar_first_child_jump'],theta['hbar_child_rooms']
    plan['held_coordinates']={r['parameter']:float(r['value']) for r in fixed.read_csv(Path(rows[0]['summary']).parent/'parameter_table.csv') if r['is_free_parameter']=='True' and r['parameter']!='hbar_child_rooms'}
    if len(plan['held_coordinates'])!=9:raise RuntimeError('Expected nine held coordinates')
    plan['preserved_three_child_floor']=j0+3*c0
    center=fixed.read_json(OUT/'profile_smoke'/prior['cases'][0]['center'])
    upper=min(.9,j0+3*(c0-.1))
    if upper<.75:raise RuntimeError('Declared extension would violate per-child lower bound')
    for jump in sorted(set((.6,.75,upper))):
        slope=max(.1,c0-(jump-j0)/3)
        if abs(jump+3*slope-(j0+3*c0))>1e-12:raise RuntimeError('Three-child floor changed')
        payload={'old_psi_child':center['old_psi_child'],'best_candidate':{
            'theta':copy.deepcopy(center['best_candidate']['theta']),
            'new_psi_child':center['best_candidate']['new_psi_child']}}
        payload['best_candidate']['theta'].update(hbar_first_child_jump=jump,hbar_child_rooms=slope)
        payload['proposal']={'status':'unevaluated conditional profile, not a calibration','fixed_jump':jump,'three_child_floor':j0+3*c0}
        planner.add_case(plan,out,payload,f'fixed_jump_{jump:.6f}_preserve_three',radius=.0025)
        plan['cases'][-1]['fixed_first_child_room_jump']=jump
    planner.finish_plan(plan,out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('smoke','grid','repeat','collect'))
    p.add_argument('--prior',type=Path);p.add_argument('--prior-sha256');p.add_argument('--allow-partial',action='store_true')
    a=p.parse_args()
    if a.action=='smoke':smoke()
    else:
        prior,status,rows=collect(a.prior,a.prior_sha256,complete=not a.allow_partial)
        if a.action=='grid':grid(prior,rows)
        elif a.action=='repeat':
            if prior['stage']!='conditional_profile_three_child_floor_preserved':raise RuntimeError('Complete profile grid first')
            repeated_plan(prior,fixed.read_json(OUT/'profile_best_so_far.json')['best'],OUT/'profile_repeats','conditional_profile_final_exact_repeats')
        else:print(json.dumps(status))
