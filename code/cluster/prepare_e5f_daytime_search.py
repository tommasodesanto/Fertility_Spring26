#!/usr/bin/env python3
"""Prepare a fresh Torch continuation; authenticate unchanged model and targets."""
import copy, hashlib, json, pathlib, time
ROOT=pathlib.Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a')
OUT=ROOT/'daytime_calibration_20260927/search'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p): return dict(path=str(p),sha256=sha(p))
def canonical(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write(p,x):
    with p.open('x') as f: json.dump(x,f,indent=2,sort_keys=True);f.write('\n')
base=ROOT/'calibration_code_integration_20260927_v2/launch_v1/contract.json'
c=read(base);old=read(OUT/'input/overnight_contract.json');obj=read(OUT/'input/overnight_objective.json');receipt=read(OUT/'input/selected_receipt.json')
assert sha(base)=='399abb6e9eab0d447dca627f94e3de4e6a8006d8920d02241fadac21f6c6ebae'
assert sha(OUT/'input/overnight_contract.json')=='3b770d8c8c22d2b0449b34a575d6353b063bc015d74ce11016dad7e22ed7ca5e'
assert receipt['status']=='verified_provisional_calibration_point'
assert receipt['source_manifest_sha256']==c['source_manifest']['sha256']==old['source_manifest']['sha256']
assert receipt['target_system_sha256']==sha(OUT/'input/overnight_objective.json')==old['objective']['sha256']
remote_obj=read(pathlib.Path(c['objective']['path']))
assert remote_obj['target_rows']==obj['target_rows']
assert remote_obj['parameter_restrictions']==obj['parameter_restrictions']
assert canonical(obj['target_rows'])==old['target_weight_fingerprint']
manifest=read(pathlib.Path(c['source_manifest']['path']))
for rel,h in manifest['files'].items(): assert sha(pathlib.Path(c['source_root'])/rel)==h
# Restore only the original ancestry lock literal in the fresh tools copy.
runner=OUT/'tools/run_e5f_utility_comparison.py';text=runner.read_text()
local_lock='6c5e14b40eba0ca63911f5a4514d3a63832259f2f36458b05c57d801eb0dadad'
original_lock='6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1'
assert text.count(local_lock)==1
runner.write_text(text.replace(local_lock,original_lock))
assert sha(runner)==c['files']['recovery_runner']['sha256']
oldtools=pathlib.Path(c['runtime_tools']);newtools=OUT/'tools'
for k,r in list(c['files'].items()):
    p=pathlib.Path(r['path'])
    if p.is_relative_to(oldtools): c['files'][k]=pin(newtools/p.relative_to(oldtools))
    else: assert sha(p)==r['sha256']
inv={'root':str(newtools),'files':{str(p.relative_to(newtools)):sha(p) for p in sorted(newtools.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo')}}
inv['file_count']=len(inv['files']);write(OUT/'runtime_inventory.json',inv)
c['files']['runtime_inventory']=pin(OUT/'runtime_inventory.json')
for rel in inv['files']: c['files']['runtime_file:'+rel]=pin(newtools/rel)
c['files']['daytime_preparer']=pin(pathlib.Path(__file__))
for name in ['overnight_contract.json','overnight_objective.json','selected_receipt.json','selected_target_fit.csv','selected_parameters.csv']: c['files']['daytime_input:'+name]=pin(OUT/'input'/name)
c.update(runtime_tools=str(newtools),initial_point=receipt['point'],execution={'kind':'slurm'},seed=2026092711,objective=pin(OUT/'input/overnight_objective.json'),objective_canonical_sha256=canonical(obj),target_weight_fingerprint=old['target_weight_fingerprint'])
for key in ['fixed','normalization','proposal_widths','identification','pending_observer_mismatches','economic_changes']:
    assert c[key]==old[key] or key in ['identification','economic_changes']
    c[key]=copy.deepcopy(old[key])
c['budget'].update(workers=24,points_per_round=24,rounds=8,absolute_end_epoch=time.time()+12600)
c['approval']={'authority':'Tommaso September27 daytime: continue search on cluster; same model and targets','production_authorized':True,'condition':'Fresh exact-loop smoke and lead cross-host numerical acceptance'}
c['status']='reviewed_smoke';c['production_blockers']=['Fresh exact-loop smoke and lead acceptance pending']
c['scope']='Bounded daytime continuation from overnight main de_0093; no economic, grid, target, weight, bound, normalization, or gate changes'
write(OUT/'contract.json',c)
write(OUT/'preparation_receipt.json',dict(contract=pin(OUT/'contract.json'),source_manifest=c['source_manifest'],target_weight_fingerprint=c['target_weight_fingerprint'],initial_point=c['initial_point'],initial_loss=receipt['loss'],workers=24,max_search_cases=192,absolute_end_epoch=c['budget']['absolute_end_epoch'],changes=['initial_point','random_seed','24 workers','8 rounds','absolute 3.5-hour cap','Torch execution and original ancestry paths'],model_source_bytes_unchanged=True,target_rows_unchanged=True,parameter_bounds_unchanged=True))
print(json.dumps(read(OUT/'preparation_receipt.json'),indent=2))
