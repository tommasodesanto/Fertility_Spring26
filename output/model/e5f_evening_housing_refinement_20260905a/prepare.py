#!/usr/bin/env python3
"""Prepare the predeclared evening plans; no model is solved here."""
import argparse
import copy
import datetime
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'code/model'), str(ROOT/'code/model/tools')]
import run_e5f_bounded_calibration_refinement as adapter
import build_e5f_bounded_refinement_plan as planner


def smoke():
    original = ROOT/'output/model/e5f_morning_refinement_20260905a'
    out = OUT/'smoke'
    out.mkdir(exist_ok=False)
    reference = ROOT/'output/model/e5f_policy_cleanup_verification_20260905a/history'
    center = original/'repeats/center_001.json'
    adapter.verify(center, 'bb8ee43280b44ad1f01acc70a2c6f470017605c22ea0d5feb10a5d540b44e5a5')
    ref = out/'reference'
    ref.mkdir()
    summary = adapter.read_json(reference/'summary.json')
    adapter.verify(reference/'summary.json', 'd7a376805b2decf2ea56be65d4ad0f11c2d18a7039befc61b697141fdaaad9c8')
    for name in ('summary.json', 'target_fit_long.csv', 'parameter_table.csv'):
        shutil.copy2(reference/name, ref/name)
    name = summary['best_candidate']['candidate']
    (ref/'cases'/name).mkdir(parents=True)
    shutil.copy2(reference/'cases'/name/'transition_path.csv', ref/'cases'/name/'transition_path.csv')
    plan = adapter.read_json(original/'smoke/plan.json')
    plan.update(stage='evening_two_exact_smokes', cases=[],
        code_bundle_sha256=adapter.VERIFIED_REPAIR_BUNDLE,
        launch_deadline_epoch=datetime.datetime(2026,9,6,0,10,tzinfo=datetime.timezone.utc).timestamp(),
        adapter_sha256=adapter.digest(adapter.__file__),
        planner_sha256=adapter.digest(planner.__file__),
        evening_preparation_sha256=adapter.digest(__file__),
        helper_sha256={name:adapter.digest(ROOT/'tmp/e5f_policy_cleanup_20260905a/code/model/tools'/name)
            for name in ('collect_e5f_transition_calibration.py', 'run_e5f_independent_numerical_audit.py')},
        input_sha256={str(p.relative_to(out)):adapter.digest(p) for p in ref.rglob('*') if p.is_file()})
    for i in range(2):
        planner.add_case(plan,out,adapter.read_json(center),f'candidate_exact_repeat_{i+1}',radius=.005)
        plan['cases'][-1].update(reference='reference/summary.json',reference_sha256=adapter.digest(ref/'summary.json'))
    planner.finish_plan(plan,out)


def coordinate(prior, status, rows, out):
    if len(rows)!=2 or rows[0]['loss']!=rows[1]['loss']:
        raise RuntimeError('Two exact smoke repeats required')
    for row in rows:
        receipt=adapter.read_json(Path(row['summary']).parent/'case_receipt.json')
        if not receipt.get('reference',{}).get('exact_twelve_row_fit'):
            raise RuntimeError('Missing exact smoke fit')
    plan=planner.fresh_plan(prior,out,'evening_coordinate_0025')
    center=adapter.read_json(rows[0]['summary'])
    for i in range(1,24):
        planner.add_case(plan,out,center,'anchor' if i==1 else f'coordinate_{(i-2)//2}_{"minus" if i%2==0 else "plus"}',
            panel_id=i,panel_size=23,panel_design='coordinate',radius=.0025)
    planner.finish_plan(plan,out)


def joint(prior,status,rows,out):
    if prior['stage']!='evening_coordinate_0025' or len(rows)!=23:
        raise RuntimeError('Complete evening coordinate panel required')
    import run_e5f_transition_calibration as calibration
    by_id={r['id']:r for r in rows}
    best=adapter.read_json(rows[0]['summary'])
    domain=best['panel_design']['domain']
    names=[d['name'] for d in domain]
    units=best['panel_design']['unit_vector']
    theta=best['best_candidate']['theta']
    excluded={'H0','hbar_child_rooms','hbar_first_child_jump'}
    direction=[]
    for j,name in enumerate(names):
        winner=min((by_id[2+2*j],by_id[3+2*j]),key=lambda r:(r['loss'],r['id']))
        improving=winner['loss']<by_id[1]['loss']-1e-8 and name not in excluded
        direction.append((-.0025 if winner['id']%2==0 else .0025) if improving else 0.)
    jump=next(d for d in domain if d['name']=='hbar_first_child_jump')['upper']
    delta=jump-theta['hbar_first_child_jump']
    patterns=[]
    for preserve in (3,2,None):
        for multiplier in (1.,.97,.85):
            changes={'hbar_first_child_jump':jump,'H0':theta['H0']*multiplier,
                'hbar_child_rooms':theta['hbar_child_rooms']-(delta/preserve if preserve else 0.)}
            patterns.append((f'jump_bound_preserve_{preserve}_H0_{multiplier}',changes,0.))
    for scale in (.5,1.,2.):
        patterns.append((f'jump_bound_preserve_3_successful_moves_{scale}',
            {'hbar_first_child_jump':jump,'hbar_child_rooms':theta['hbar_child_rooms']-delta/3},scale))
    plan=planner.fresh_plan(prior,out,'evening_joint_housing_endpoint')
    plan.update(starting_summary_sha256=rows[0]['summary_sha256'],coordinate_direction=direction,
        rule='Original twelve-row objective; predeclared floor endpoint and supply patterns; no Jacobian inverse')
    seen={tuple(units)}
    for label,changes,scale in patterns:
        unit=[min(1.,max(0.,u+scale*d)) for u,d in zip(units,direction)]
        for name,value in changes.items():
            j=names.index(name);spec=domain[j]
            if not spec['lower']<=value<=spec['upper']:
                raise RuntimeError(f'Predeclared housing pattern leaves bounds: {name}')
            unit[j]=calibration.inverse_unit(value,spec['lower'],spec['upper'],spec['transform'])
        if tuple(unit) in seen:
            continue
        seen.add(tuple(unit))
        payload={'old_psi_child':best['old_psi_child'],'best_candidate':{
            'theta':copy.deepcopy(theta),'new_psi_child':best['best_candidate']['new_psi_child'],
            'old_psi_child':best['old_psi_child']}}
        for u,spec in zip(unit,domain):
            value=calibration.transform_unit(u,spec['lower'],spec['upper'],spec['transform'])
            name=spec['name']
            if name=='beta_annual':payload['best_candidate']['theta']['beta']=value**4
            elif name=='psi_child_change_2023':payload['best_candidate']['new_psi_child']=best['old_psi_child']+value
            else:payload['best_candidate']['theta'][name]=value
        payload['proposal']={'label':label,'unit_vector':unit,'status':'unevaluated, no fitted moments or loss'}
        planner.add_case(plan,out,payload,label,radius=.0025)
    if len(plan['cases'])>12:raise RuntimeError('Joint budget exceeded')
    planner.finish_plan(plan,out)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('smoke','coordinate','joint','repeat'))
    parser.add_argument('--prior',type=Path)
    parser.add_argument('--prior-sha256')
    args=parser.parse_args()
    if args.action=='smoke':smoke()
    else:
        prior,status,rows=planner.collect(args.prior,args.prior_sha256)
        adapter.verify(__file__,prior['evening_preparation_sha256'])
        stage={'coordinate':'round1','joint':'round2','repeat':'repeats'}[args.action]
        {'coordinate':coordinate,'joint':joint,'repeat':planner.repeats}[args.action](prior,status,rows,OUT/stage)
