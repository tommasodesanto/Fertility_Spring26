"""Torch-side pinned manifests for the September13 rebated finite histories."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
INITIAL=Path('/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8')
BATCH=ROOT/'batches/final_night_20260913'

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--initial-summary',type=Path,required=True)
    parser.add_argument('--label',default='joint')
    args=parser.parse_args()
    if not args.label.replace('_','').isalnum():raise ValueError('Simple manifest label required')
    source=BATCH/'history_source';pairs=[]
    for first in sorted((INITIAL/'code/model/intergen_eqscale_seq_optimized').rglob('*.py')):
        second=ROOT/first.relative_to(INITIAL)
        digest=sha(first)
        if digest!=sha(second):raise ValueError('Initial/history economic kernel mismatch: '+str(first))
        pairs.append(dict(initial=str(first),history=str(second),sha256=digest))
    equivalence=BATCH/'kernel_equivalence.json'
    save(equivalence,dict(status='equal',pairs=pairs))
    prior=json.loads((ROOT/'batches/night_surprises_20260912/plan.json').read_text())
    prior['maximum_trials_per_window']=6
    prior['launch_refinement']='Six trials per window; explicit rebated finite boundaries; original target and source pins retained'
    plan=BATCH/'history_plan.json';save(plan,prior)
    pins={str(f):sha(f) for f in source.glob('*.py')}
    pins.update({str(plan):sha(plan),str(equivalence):sha(equivalence),prior['empirical_blocks']:sha(prior['empirical_blocks'])})
    forecast=BATCH/f'history_manifest_{args.label}.json'
    save(forecast,dict(prior_plan=str(plan),initial_summary=str(args.initial_summary),
        initial_source_root=str(INITIAL),kernel_equivalence=str(equivalence),file_sha256=pins,
        policy_reserve_seconds=10800,seed_step=-.01,
        root_controls=dict(transfer_bounds=[1e-10,10.]),
        disclosure='A historical observed head-age conditioning; dated equal rebate and PAYGO; finite boundary actual carried state, future constant conditions provisional'))
    stages=[]
    for case,label in [('A0','A0'),('A+','Aplus')]:
        for count in (6,24,100):
            name=f'{label}_{count}'
            stages.append(dict(name=name,seconds=43200,cpus=1,
                command=['python','-B',str(source/'run_e5f_final_rebated_history.py'),
                    '--manifest',str(forecast),'--case',case,'--count',str(count),
                    '--output',str(BATCH/f'histories_{args.label}'/name),'--seconds','42600']))
    submission=BATCH/f'history_array_manifest_{args.label}.json'
    save(submission,dict(runroot=str(BATCH/f'histories_{args.label}'),horizon_hours=12,max_workers=6,memory='32G',
        source_pins={**pins,str(forecast):sha(forecast)},stages=stages,
        environment=dict(PYTHONPATH=f'{source}:{ROOT}/code/model/tools:{ROOT}/code/model',
            NUMBA_CACHE_DIR=str(ROOT/'output/cache/numba'))))
    print(json.dumps(dict(manifest=str(submission),forecast_manifest=str(forecast),kernel_files=len(pairs),
        tracks=len(stages),B_status='formation and migration allocation outstanding; not submitted as calibrated histories')))

if __name__=='__main__':main()
