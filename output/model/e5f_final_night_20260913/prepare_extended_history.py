"""Torch-only, pinned 24-to-100 numerical guesses; no imported household state."""
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time

B=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913')
old=B/'history_source_seeded_v1'; new=B/'history_source_extended_seed_v1'
auto=B/'history_source_auto_v1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2);f.write('\n')

for p in old.glob('*.py'):
    if p.name not in ('run_e5f_final_rebated_history.py','test_run_e5f_final_rebated_history.py'):
        shutil.copyfile(p,new/p.name)
assert sorted(p.name for p in new.glob('*.py') if sha(p)!=sha(old/p.name))==[
    'run_e5f_final_rebated_history.py','test_run_e5f_final_rebated_history.py']
base_path=B/'history_manifest_corrected_seeded_v1.json'
base=read(base_path); manifest=copy.deepcopy(base)
pins={str(new/Path(p).name) if Path(p).parent==old else p:h for p,h in base['file_sha256'].items()}
for p in new.glob('*.py'):pins[str(p)]=sha(p)
for key,case in (('A0','A0'),('Aplus','A+')):
    root=B/'histories_refit'/f'{key}_24/window_2007/trial_00/root_receipt.json'
    r=read(root);prices=r['best']['prices']
    assert r['case']==case and r['count']==24 and r['start_year']==2007
    assert r['best']['mapping_valid'] is True and len(prices)==75
    coordinates=[]
    for start in (0,25,50):
        block=prices[start:start+25];coordinates.extend(block+[block[-1]]*76)
    assert len(coordinates)==303 and all(math.isfinite(x) and x>0 for x in coordinates)
    path=B/'numerical_coordinate_seeds'/f'{key}_100_extended24.json'
    save(path,dict(label='numerical_guess_only',case=case,count=100,start_year=2007,
        selection='best',coordinates=coordinates,source_root_receipt=dict(path=str(root),sha256=sha(root)),
        source_count=24,rule='hold_last'))
    manifest['initial_coordinate_seeds'][f'{key}_100']=dict(path=str(path),sha256=sha(path))
    pins[str(path)]=sha(path);pins[str(root)]=sha(root)
manifest['file_sha256']=dict(sorted(pins.items()))
comparison=copy.deepcopy(manifest)
comparison['file_sha256']=base['file_sha256'];comparison['initial_coordinate_seeds']=base['initial_coordinate_seeds']
assert comparison==base
destination=B/'history_manifest_corrected_extended_seed_v1.json';save(destination,manifest)
sys.path.insert(0,str(new))
from run_e5f_final_rebated_history import initial_coordinate_seed, verify_pins
import numpy as np
verify_pins(pins)
preflight={case:initial_coordinate_seed(manifest,case,100,np.ones(303))[1] for case in ('A0','A+')}
remaining=math.floor(1789322400-time.time());assert 10860<remaining<=43200
parent=B/'history_array_manifest_corrected_cpu32_100.json';array=read(parent)
runroot=B/'histories_corrected_extended_100';array['runroot']=str(runroot)
array['horizon_hours']=remaining/3600
array['environment']['PYTHONPATH']=array['environment']['PYTHONPATH'].replace(str(auto),str(new))
array_pins={p:h for p,h in array['source_pins'].items() if Path(p).parent!=auto and p!=str(B/'history_manifest_corrected_auto_v2.json')}
array_pins.update(pins);array_pins[str(destination)]=sha(destination)
array['source_pins']=dict(sorted(array_pins.items()))
for stage in array['stages']:
    argv=stage['command'];stage['seconds']=remaining
    for flag,value in [('--manifest',str(destination)),('--output',str(runroot/stage['name'])),('--seconds',str(remaining-60))]:
        assert argv.count(flag)==1;argv[argv.index(flag)+1]=value
path=B/'history_array_manifest_corrected_extended_100.json';save(path,array)
receipt=dict(status='prepared_not_submitted',parent=str(parent),parent_sha256=sha(parent),
    manifest=str(destination),manifest_sha256=sha(destination),array_manifest=str(path),array_sha256=sha(path),
    driver_sha256=sha(new/'run_e5f_final_rebated_history.py'),preflight=preflight,
    numerical_guess_only=True,household_state_reuse=False,jacobian_reuse_from_seed=False,
    structural_parameters_changed=False,targets_or_gates_changed=False,
    sizing='Two independent100period histories; four windows, at most six preference trials perwindow and24root mappings pertrial plus retained bounded alternative. Recent100period nonconstant mapping exceeds10minutes; finishing allwindows before18UTC is not assured. Existing case checkpoints, latest/best summaries, firstnative mapping gate,30minute health threshold and twohour policy reserve retained.')
save(B/'extended_history_preparation.json',receipt)
print(json.dumps({k:v for k,v in receipt.items() if k!='preflight'}))
