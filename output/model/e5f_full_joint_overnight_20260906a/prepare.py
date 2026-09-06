#!/usr/bin/env python3
"""Prepare the immutable overnight contract and remote-path seed, without solves."""
from pathlib import Path
import copy
import datetime
import hashlib
import json
import math
import shutil
import sys
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
SNAP=ROOT/'tmp/e5f_full_joint_overnight_20260906a'
REMOTE=Path('/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906a')
REL=OUT.relative_to(ROOT)
sys.path.insert(0,str(ROOT/'code/model/tools'))
import run_e5f_joint_overnight_case as adapter

def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def inverse(value,spec):
    low,high,kind=spec['lower'],spec['upper'],spec['transform']
    if not low<=value<=high:raise ValueError((value,spec))
    if kind=='log':return math.log(value/low)/math.log(high/low)
    if kind=='discount':return 1-math.sqrt(max(0,1-(value-low)/(high-low)))
    if kind=='asinh':return (math.asinh(value)-math.asinh(low))/(math.asinh(high)-math.asinh(low))
    if kind=='softzero':return math.sqrt((value-low)/(high-low))
    raise ValueError(kind)

def main():
    if (OUT/'contract.json').exists():raise RuntimeError('Immutable contract already exists')
    prior=ROOT/'output/model/e5f_evening_housing_refinement_20260905a'
    inputs=OUT/'inputs';inputs.mkdir()
    seed=json.loads((prior/'empirical_room_gap/grid/center_004.json').read_text())
    delta=seed['best_candidate']['new_psi_child']-seed['old_psi_child']
    values=[]
    for spec in adapter.SEARCH_DOMAIN:
        name=spec['name']
        value=delta if name=='psi_child_change_2023' else seed['best_candidate']['theta']['beta']**.25 if name=='beta_annual' else seed['best_candidate']['theta'][name]
        values.append(inverse(value,spec))
    seed['panel_design']={'unit_vector':values,'domain':adapter.SEARCH_DOMAIN}
    save(inputs/'seed_center.json',seed)
    ref=prior/'empirical_room_gap/grid/task_004';dest=inputs/'reference';dest.mkdir()
    summary=json.loads((ref/'summary.json').read_text());case=summary['best_candidate']['candidate']
    for name in ('summary.json','target_fit_long.csv','parameter_table.csv'):
        shutil.copy2(ref/name,dest/name)
    (dest/'cases'/case).mkdir(parents=True)
    shutil.copy2(ref/'cases'/case/'transition_path.csv',dest/'cases'/case/'transition_path.csv')
    source=json.loads((OUT/'source_contract.json').read_text())
    base=json.loads((prior/'smoke/plan.json').read_text())
    base={k:base[k] for k in ('source','source_sha256','target_set','target_fingerprint','helper_sha256')}
    base.update(schema='e5f_joint_overnight_v1',cases=[],input_sha256={},production_promoted=False,
        code_bundle_sha256=source['bundle_sha256'],adapter_sha256=adapter.digest(adapter.__file__),
        search_domain=adapter.SEARCH_DOMAIN,first_child_jump_upper=2.,
        planner_sha256=adapter.digest(ROOT/'code/model/tools/build_e5f_bounded_refinement_plan.py'),
        launch_deadline_epoch=datetime.datetime(2026,9,6,11,tzinfo=datetime.timezone.utc).timestamp())
    c=dict(schema='e5f_joint_overnight_contract_v1',base_plan=base,search_domain=adapter.SEARCH_DOMAIN,
        output_root=str(REMOTE/REL),seed_center=str(REMOTE/REL/'inputs/seed_center.json'),
        seed_reference=str(REMOTE/REL/'inputs/reference/summary.json'),
        seed_center_sha256=adapter.digest(inputs/'seed_center.json'),seed_reference_sha256=adapter.digest(dest/'summary.json'),
        controller_sha256=adapter.digest(ROOT/'code/model/tools/run_e5f_joint_overnight_search.py'),
        case_adapter_sha256=adapter.digest(adapter.__file__),planner_sha256=base['planner_sha256'],
        source_sha256=base['source_sha256'],target_fingerprint=base['target_fingerprint'],code_bundle_sha256=source['bundle_sha256'],
        seed=2026090601,max_workers=12,case_timeout_seconds=1800,max_histories=360,
        max_search_seconds=27000,max_total_seconds=30600,max_generations=8,population_size=32,
        initial_radius=.01,polish_rounds=2,polish_radius=.0025,
        absolute_finish_epoch=datetime.datetime(2026,9,6,12,tzinfo=datetime.timezone.utc).timestamp(),
        smoke_dir=str(REMOTE/REL/'smoke'),production_promoted=False)
    save(OUT/'contract.json',c)
    for name in ('run_e5f_joint_overnight_search.py','run_e5f_joint_overnight_case.py',
                 'test_run_e5f_joint_overnight_search.py','build_e5f_bounded_refinement_plan.py'):
        shutil.copy2(ROOT/'code/model/tools'/name,SNAP/'code/model/tools'/name)
    print(json.dumps({'contract_sha256':adapter.digest(OUT/'contract.json'),'bundle':source['bundle_sha256'],'estimated_parameters':11,'targets':12}))

if __name__=='__main__':main()
