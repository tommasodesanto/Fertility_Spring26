"""Prepare one immutable, deadline-bounded smoke or historical timing case."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import run_e5f_matched_pf_preference_pilot as pilot

PARENT = Path('/scratch/td2248/projects/Fertility_Spring26_matched_pf_root_h100_20260910a/output/historical_root_h100_continue_01/sequential')
PINS = {'contract.json':'442dd06ebe369ca1e282963a6926effede1775dcca0fceb135b070e50d26060b',
    'root_history.json':'7d7b6d7f78eaabf8aa5a87df64febfaa14d502fc047e0894e2803270dc772020',
    'summary.json':'915eeeee6dbe24b3fc45bd1517a647c9188d765fa8c3886f45e28ab447b80284'}
DEADLINE = 1789086600  # 2026-09-11 00:30 UTC; checked before every case starts.


def check_completed_smokes(root):
    for i,shape in enumerate(pilot.SHAPES):
        out=root/'output'/'smoke'/f'case_{i}'
        s=json.loads((out/'summary.json').read_text())
        if (s['status']!='passed_conditional_preference_pilot' or s['shape']!=shape
                or s['scope']!='six_date_plumbing_smoke_only' or s['evaluations']!=1):
            raise ValueError('All three exact case smokes must pass before a main run')
        e=json.loads((out/'evaluation_001'/'summary.json').read_text())
        if not e['mapping_valid'] or not all(g['passed'] for g in e['gates'].values()):
            raise ValueError('Smoke numerical or accounting gate failed')
        for name,expected in e['artifact_sha256'].items():
            pilot.primitive.verify(out/'evaluation_001'/name,expected)
        with (out/'evaluation_001'/'transition_path.csv').open() as stream:
            years=[int(row['calendar_year']) for row in csv.DictReader(stream)]
        if years!=[2007,2011,2015,2019,2023,2027]:
            raise ValueError('Smoke dates incomplete')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase',required=True,choices=('smoke','main'))
    p.add_argument('--case',required=True,type=int,choices=(0,1,2))
    p.add_argument('--validate-only',action='store_true')
    args=p.parse_args()
    root=Path(__file__).resolve().parents[3]
    remaining=int(DEADLINE-time.time())-60
    if remaining<120:
        raise TimeoutError('Author review-window budget is exhausted; no solve started')
    for name,pin in PINS.items():pilot.primitive.verify(PARENT/name,pin)
    parent=json.loads((PARENT/'contract.json').read_text())
    c=copy.deepcopy(parent)
    for name in list(c):
        if name.startswith('restart_'):del c[name]
    for field in ('contract_sha256','numerical_controls'):
        c.pop(field,None)
    shape=pilot.SHAPES[args.case]
    smoke=args.phase=='smoke'
    mode='conditional' if smoke else ('replay' if shape==0 else 'root')
    c.update(experiment='fixed_endpoint_preference_shape_pilot',arm='sequential',
        shape_coefficient=shape,pilot_mode=mode,probe_coordinate=-1,path_date_count=6 if smoke else 100,
        preconditioner='diagonal' if smoke else 'verified_parent_broyden',
        maximum_path_evaluations=3 if mode=='root' else 1,
        seconds=min(1200 if smoke else 10800,remaining),absolute_deadline_unix=DEADLINE,
        scope='six_date_plumbing_smoke_only' if smoke else 'fixed_endpoint_historical_preference_timing_pilot',
        save_initial_2023_state=False)
    for field,name in [('parent_contract','contract.json'),('parent_history','root_history.json'),('parent_summary','summary.json')]:
        c[field]=str(PARENT/name);c[field+'_sha256']=PINS[name]
    driver=pilot.DRIVER
    for name,expected in parent['source_sha256'].items():
        if name!=driver:pilot.primitive.verify(root/name,expected)
    c['source_sha256'][driver]=pilot.primitive.digest(root/driver)
    for filename in ('run_e5f_matched_pf_preference_pilot.py','test_run_e5f_matched_pf_preference_pilot.py',
            'prepare_e5f_matched_pf_preference_case.py','e5f_matched_pf_birth_path.py','test_e5f_matched_pf_birth_path.py'):
        path='code/model/tools/'+filename
        c['source_sha256'][path]=pilot.primitive.digest(root/path)
    c['reviewed_preference_driver_change']=dict(path=driver,
        from_sha256=parent['source_sha256'][driver],to_sha256=c['source_sha256'][driver],scope=pilot.CHANGE_SCOPE)
    pilot.load_parent(c,'sequential')
    if not smoke and not args.validate_only:check_completed_smokes(root)
    if args.validate_only:
        print(json.dumps(dict(status='validated_without_solve',phase=args.phase,case=args.case,
            shape=shape,mode=mode,seconds=c['seconds'],source_files=len(c['source_sha256']))))
        return
    folder=root/'contracts'/args.phase
    folder.mkdir(parents=True,exist_ok=True)
    contract=folder/f'case_{args.case}.json'
    if contract.exists():raise FileExistsError('Refusing duplicate case contract')
    contract.write_text(json.dumps(c,indent=2)+'\n')
    pin=pilot.primitive.digest(contract)
    output=root/'output'/args.phase/f'case_{args.case}'
    pilot.joined.load_smoke_contract(contract,pin,'sequential',maximum_seconds=10800)
    subprocess.run([sys.executable,str(root/'code/model/tools/run_e5f_matched_pf_preference_pilot.py'),
        '--contract',str(contract),'--contract-sha256',pin,'--arm','sequential',
        '--shape',str(shape),'--mode',mode,'--output',str(output)],check=True)


if __name__=='__main__':main()
