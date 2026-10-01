"""Bounded, read-only reduction of compact checkpoint JSONs; never imports model code."""
import argparse, json, statistics, time
from collections import Counter
from pathlib import Path

def collect(root, remote=False):
    rows=[]
    for folder in sorted(root.glob('chain_*')):
        search=folder/('search' if remote else 'results')
        def read(name):
            p=search/name
            return json.loads(p.read_text()) if p.exists() else None
        cases=read('cases.json') or []; passed=[r for r in cases if r.get('status')=='passed']; contract=read('input_contract.json') or {}; latest=read('latest.json') or {}
        def best(key):return min(passed,key=lambda r:r.get(key,r.get('loss',float('inf')))) if passed else None
        base=best('base_loss');weighted=best('objective');times=[r['seconds'] for r in passed if isinstance(r.get('seconds'),(int,float))]
        active=search/(latest.get('label') or '')/'latest.json';native=json.loads(active.read_text()) if active.is_file() else None
        rows.append(dict(chain=folder.name,profile=contract.get('weight_profile'),completed_cases=len(cases),admissible_computed_ge=len(passed),statuses=dict(Counter(r.get('status') for r in cases)),objective_calls=latest.get('objective_calls'),native_lifecycle_total=sum(r.get('lifecycle_solves',0) for r in cases),case_seconds=dict(count=len(times),minimum=min(times),median=statistics.median(times),maximum=max(times)) if times else None,best_base=base,best_weighted=weighted,latest=latest,native_latest=native,failure=read('failure.json'),completed=read('completed.json'),postcheck=(json.loads((folder/'postcheck/completed.json').read_text()) if (folder/'postcheck/completed.json').exists() else None)))
    return dict(checked_epoch=time.time(),root=str(root),remote=remote,rows=rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--remote',action='store_true');a=p.parse_args();print(json.dumps(collect(a.root,a.remote)))
