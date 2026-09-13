"""Torch-only numerical warm-start manifests; no household states are imported."""
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import time

B=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913')
old=B/'history_source_auto_v1';new=B/'history_source_seeded_v1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2);f.write('\n')
assert sha(new/'run_e5f_final_rebated_history.py')=='10615e4dc66c50d3ac7708a307dda6a2051134b47b99ad8e452ce7066d8b053f'
for p in old.glob('*.py'):
    if p.name!='run_e5f_final_rebated_history.py':shutil.copyfile(p,new/p.name)
assert sorted(p.name for p in new.glob('*.py') if sha(p)!=sha(old/p.name))==['run_e5f_final_rebated_history.py']
seeds=B/'numerical_coordinate_seeds'
for key in ('A0_24','Aplus_24'):
    root=B/'histories_refit'/key/'window_2007/trial_00/root_receipt.json'
    r=read(root);coords=r['best']['prices']
    assert r['count']==24 and r['start_year']==2007 and r['best']['mapping_valid'] and len(coords)==75
    assert all(math.isfinite(x) and x>0 for x in coords)
    save(seeds/(key+'.json'),dict(label='numerical_guess_only',case=r['case'],count=24,start_year=2007,
        selection='best',coordinates=coords,source_root_receipt=dict(path=str(root),sha256=sha(root))))
base_manifest=B/'history_manifest_corrected_auto_v2.json';destination=B/'history_manifest_corrected_seeded_v1.json'
base=read(base_manifest);manifest=copy.deepcopy(base)
pins={str(new/Path(p).name) if Path(p).parent==old else p:h for p,h in base['file_sha256'].items()}
for p in new.glob('*.py'):pins[str(p)]=sha(p)
profiles={}
for p in sorted(seeds.glob('*.json')):
    d=read(p);origin=d['source_root_receipt'];assert sha(origin['path'])==origin['sha256']
    profiles[p.stem]=dict(path=str(p),sha256=sha(p));pins[str(p)]=sha(p);pins[origin['path']]=origin['sha256']
manifest['file_sha256']=dict(sorted(pins.items()));manifest['initial_coordinate_seeds']=profiles
comparison=copy.deepcopy(manifest);comparison['file_sha256']=base['file_sha256'];comparison.pop('initial_coordinate_seeds')
assert comparison==base
save(destination,manifest)
remaining=math.floor(1789322400-time.time());assert 10860<remaining<=43200
receipts={}
for count in (6,24,100):
    array=read(B/f'history_array_manifest_corrected_auto_v2_{count}.json')
    runroot=B/f'histories_corrected_seeded_v1_{count}'
    array['runroot']=str(runroot);array['horizon_hours']=remaining/3600
    array['environment']['PYTHONPATH']=array['environment']['PYTHONPATH'].replace(str(old),str(new))
    source_pins={p:h for p,h in array['source_pins'].items() if Path(p).parent!=old and p!=str(base_manifest)}
    source_pins.update(pins);source_pins[str(destination)]=sha(destination)
    array['source_pins']=dict(sorted(source_pins.items()))
    for stage in array['stages']:
        argv=stage['command'];stage['seconds']=remaining
        for flag,value in [('--manifest',str(destination)),('--output',str(runroot/stage['name'])),('--seconds',str(remaining-60))]:
            assert argv.count(flag)==1;argv[argv.index(flag)+1]=value
    path=B/f'history_array_manifest_corrected_seeded_v1_{count}.json';save(path,array)
    receipts[str(count)]=dict(path=str(path),sha256=sha(path),runroot=str(runroot))
save(B/'seeded_history_preparation.json',dict(status='prepared_not_submitted',manifest=str(destination),manifest_sha256=sha(destination),arrays=receipts,
    numerical_guess_only=True,household_state_reuse=False,structural_parameters_changed=False,targets_or_gates_changed=False))
print(json.dumps(dict(status='prepared',seed_profiles=list(profiles),arrays=receipts)))
