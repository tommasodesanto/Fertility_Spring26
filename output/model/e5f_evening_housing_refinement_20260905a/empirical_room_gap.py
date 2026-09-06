#!/usr/bin/env python3
"""One bounded follow-up: adjust the slope for the observed parent-room gap.

No scientific source or target changes. Reuse the twice-smoked fixed-jump loop.
"""
import argparse,copy,datetime,json,math,sys,time
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];RUN=OUT/'empirical_room_gap'
sys.path[:0]=[str(OUT),str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import fixed_jump_adapter as adapter
import fixed_profile as profile
import build_e5f_bounded_refinement_plan as planner
profile.OUT=RUN
MOMENT='prime30_55_parent_3plus_minus_1to2_mean_rooms'
DEADLINE=datetime.datetime(2026,9,6,0,55,tzinfo=datetime.timezone.utc).timestamp()

def new_plan(prior,out,stage):
    p=planner.fresh_plan(prior,out,stage)
    p.update(schema='e5f_fixed_first_child_profile_v1',scope='conditional_profile_not_calibration',
        production_promoted=False,adapter_sha256=adapter.digest(adapter.__file__),
        profile_preparation_sha256=adapter.digest(profile.__file__),
        empirical_gap_preparation_sha256=adapter.digest(__file__),launch_deadline_epoch=DEADLINE)
    p.pop('held_coordinates',None);p.pop('preserved_three_child_floor',None)
    return p

def collect(path,sha,complete=True):
    plan=adapter.load_plan(path,sha);adapter.verify(__file__,plan['empirical_gap_preparation_sha256'])
    if len(plan['held_nine'])!=9:raise RuntimeError('Nine held coordinates required')
    for case in plan['cases']:
        out=path.parent/case['output'];receipt=out/'case_receipt.json'
        if not receipt.exists():continue
        r=adapter.read_json(receipt);adapter.verify(out/'parameter_table.csv',r['artifact_sha256']['parameter_table.csv'])
        values={v['parameter']:float(v['value']) for v in adapter.read_csv(out/'parameter_table.csv')}
        for key,v in plan['held_nine'].items():
            if not math.isclose(values[key],v,rel_tol=0,abs_tol=2e-12):raise RuntimeError(f'Held coordinate changed: {key}')
        adapter.verify(path.parent/case['center'],case['center_sha256'])
        center=adapter.read_json(path.parent/case['center'])['best_candidate']['theta']
        for key in ('hbar_child_rooms','hbar_first_child_jump'):
            if not math.isclose(values[key],center[key],rel_tol=0,abs_tol=2e-12):raise RuntimeError(f'Proposed housing coordinate changed: {key}')
    return profile.collect(path,sha,complete)

def prepare():
    if time.time()>datetime.datetime(2026,9,6,0,35,tzinfo=datetime.timezone.utc).timestamp():raise RuntimeError('New bounded start window closed')
    stop=adapter.read_json(OUT/'profile_stop.json')
    if stop['status']!='complete_no_improving_extended_point':raise RuntimeError('Previous diagnostic must be reconciled first')
    prior=adapter.read_json(OUT/'profile_smoke/plan.json')
    profile.collect(OUT/'profile_smoke/plan.json',adapter.digest(OUT/'profile_smoke/plan.json'))
    for i in (1,2):
        r=adapter.read_json(OUT/f'profile_smoke/task_{i:03}/case_receipt.json')['reference']
        if not (r['exact_twelve_row_fit'] and r['exact_parameters'] and r['exact_numeric_history_entries']==253):raise RuntimeError('Exact fixed-option smoke missing')
    anchor=OUT/'round2/task_010';s=adapter.read_json(anchor/'summary.json')
    adapter.verify(anchor/'summary.json','e7e85b5183392e8db9e3d8c7429f9286fbd385ed1d5ec647d450e7e6e12912bf')
    fit={r['moment']:r for r in adapter.read_csv(anchor/'target_fit_long.csv')}
    par={r['parameter']:float(r['value']) for r in adapter.read_csv(anchor/'parameter_table.csv')}
    domain=[r['name'] for r in adapter.read_json(OUT/'round1/task_001/summary.json')['panel_design']['domain']]
    D={};inputs={}
    for key in ('hbar_first_child_jump','hbar_child_rooms'):
        k=domain.index(key);points=[]
        for i in (2+2*k,3+2*k):
            directory=OUT/f'round1/task_{i:03}';r=adapter.read_json(directory/'case_receipt.json')
            for name in ('target_fit_long.csv','parameter_table.csv'):
                adapter.verify(directory/name,r['artifact_sha256'][name]);inputs[str((directory/name).relative_to(OUT))]=r['artifact_sha256'][name]
            f={r['moment']:float(r['model']) for r in adapter.read_csv(directory/'target_fit_long.csv')}
            t={r['parameter']:float(r['value']) for r in adapter.read_csv(directory/'parameter_table.csv')}
            points.append((f[MOMENT],t[key]))
        D[key]=(points[1][0]-points[0][0])/(points[1][1]-points[0][1])
    out=RUN/'grid';plan=new_plan(prior,out,'empirical_large_family_gap_four_points')
    plan.update(held_nine={k:v for k,v in par.items() if k in ('H0','beta_annual','chi','first_birth_fixed_cost','kappa_fert','kappa_fert_continuation','psi_child_change_2023','theta0','theta1')},
        anchor_loss=s['best_candidate']['transition_loss'],physical_derivatives=D,derivative_source_sha256=inputs,
        identifying_row=MOMENT,rule='Local slope adjustment toward the unchanged observed parent-room gap; evaluate all twelve actual moments. Linear predictions are not model results.')
    center=adapter.read_json(OUT/'profile_smoke/center_001.json')
    gap=float(fit[MOMENT]['target'])-float(fit[MOMENT]['model'])
    for jump in (.525,.55,.575,.6):
        slope=par['hbar_child_rooms']+(gap-D['hbar_first_child_jump']*(jump-.5))/D['hbar_child_rooms']
        if not .1<=slope<=1.8:raise RuntimeError('Slope violates unchanged bound')
        payload={'old_psi_child':center['old_psi_child'],'best_candidate':{'theta':copy.deepcopy(center['best_candidate']['theta']),'new_psi_child':center['best_candidate']['new_psi_child']}}
        payload['best_candidate']['theta'].update(hbar_first_child_jump=jump,hbar_child_rooms=slope)
        planner.add_case(plan,out,payload,f'jump_{jump}_observed_room_gap',radius=.0025)
        plan['cases'][-1]['fixed_first_child_room_jump']=jump
    planner.finish_plan(plan,out)

def repeat(prior,rows):
    if len(rows)!=4:raise RuntimeError('All four new points required')
    best=min(rows,key=lambda r:r['loss'])
    if best['loss']>=prior['anchor_loss']:
        adapter.write_json(RUN/'stop.json',dict(status='complete_no_improvement',best=best,production_promoted=False));return
    out=RUN/'repeats';plan=new_plan(prior,out,'empirical_room_gap_exact_repeats')
    case=next(c for c in prior['cases'] if c['id']==best['id'])
    center=adapter.read_json(RUN/'grid'/case['center']);reference=profile.copy_reference(Path(best['summary']),out)
    plan['input_sha256']={str(p.relative_to(out)):adapter.digest(p) for p in reference.parent.rglob('*') if p.is_file()}
    for i in range(2):
        planner.add_case(plan,out,center,f'exact_empirical_gap_repeat_{i+1}',radius=.0025)
        plan['cases'][-1].update(fixed_first_child_room_jump=case['fixed_first_child_room_jump'],reference=str(reference.relative_to(out)),reference_sha256=adapter.digest(reference))
    planner.finish_plan(plan,out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','collect','repeat'))
    p.add_argument('--plan',type=Path);p.add_argument('--sha256');p.add_argument('--allow-partial',action='store_true');a=p.parse_args()
    if a.action=='prepare':prepare()
    else:
        plan,status,rows=collect(a.plan,a.sha256,not a.allow_partial)
        if a.action=='repeat':repeat(plan,rows)
        else:print(json.dumps(status))
