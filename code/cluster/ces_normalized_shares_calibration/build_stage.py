"""Build immutable CES attempt 3; authenticate every frozen reporting pin first."""
import gzip, hashlib, io, json, shutil, subprocess, tarfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]; HERE=Path(__file__).resolve().parent
DEPLOY=ROOT/'output/model/experiments/ces_normalized_shares/overnight_v1/deployment'
OUT=DEPLOY/'attempt3'
PARENT=ROOT/'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/deployment/attempt2/stage.tar.gz'
PARENT_SHA='f7a8fec4ff370fd3690c0d0068ca595b75a17dd8aaac3bd47f6009ef73ecd68b'
PLAN='output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json'
ANCHOR='output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
OVERLAY='output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime'
MANIFEST='output/model/fertility_identification_20260928/fixed_reference_manifest.json'
REMOTE='/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3'

def sha(x): return hashlib.sha256(x if isinstance(x,bytes) else Path(x).read_bytes()).hexdigest()
def rel(path):
 p=Path(path).resolve()
 if not p.is_relative_to(ROOT): raise SystemExit('declared dependency outside project root: '+str(path))
 return str(p.relative_to(ROOT))
def add_tree(source,relative):
 for p in (ROOT/relative).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.nbc','.nbi'):
   source[str(p.relative_to(ROOT))]=p.read_bytes()
def historical_bytes(relative,digest):
 """Read only a named retained archive object; never restore the deleted path."""
 archived='calibration_archive/model_legacy_20261003/'+relative.removeprefix('code/model/')
 try: data=subprocess.check_output(['git','show','HEAD:'+archived],cwd=ROOT)
 except subprocess.CalledProcessError: raise SystemExit('authenticated bytes unavailable: '+relative)
 if sha(data)!=digest: raise SystemExit('retained archive hash mismatch: '+relative)
 return data
def exact_bytes(relative,digest):
 p=ROOT/relative
 if p.is_file() and sha(p)==digest: return p.read_bytes(), 'checkout'
 frozen=ROOT/OVERLAY/'frozen_sources'/Path(relative).name
 if frozen.is_file() and sha(frozen)==digest: return frozen.read_bytes(), 'frozen_overlay'
 return historical_bytes(relative,digest), 'retained_archive'
def pin_dependencies(manifest):
 """Follow only the authentication paths used by the frozen evaluator/runtime."""
 pins={}; source_maps={}
 def add_pin(pin):
  if not isinstance(pin,dict) or not {'path','sha256'}<=set(pin): return
  p=rel(pin['path']); d=pin['sha256']
  if pins.setdefault(p,d)!=d: raise SystemExit('conflicting declared digest: '+p)
 for key in ('contract','objective','source_manifest','source_contract','native_ancestry_contract'): add_pin(manifest[key])
 contract_rel=rel(manifest['contract']['path'])
 contract,_=exact_bytes(contract_rel,pins[contract_rel])
 contract=json.loads(contract)
 # authenticate() checks this exact controller list.  Its source-manifest and
 # ancestry pointers are separately verified by runtime.setup.
 for pin in contract['files'].values(): add_pin(pin)
 for key in ('source_manifest','native_ancestry_contract','source_contract'): add_pin(contract[key])
 # native.setup authenticates the ancestry contract's complete file list,
 # base contract and objective before driver.verify and runner.pair_runtime.
 ancestry_rel=rel(manifest['native_ancestry_contract']['path'])
 ancestry_data,_=exact_bytes(ancestry_rel,pins[ancestry_rel]); ancestry=json.loads(ancestry_data)
 for pin in list(ancestry['files'].values())+[ancestry['base_contract'],ancestry['objective'],ancestry['source_manifest']]: add_pin(pin)
 base_rel=rel(ancestry['base_contract']['path'])
 base_data,_=exact_bytes(base_rel,pins[base_rel]); base=json.loads(base_data)
 pair_root=Path(base['reference_root'])
 # pair_runtime pins this exact portable parent lock in the authenticated runner.
 runner_rel=rel(ancestry['files']['recovery_runner']['path'])
 runner_data,_=exact_bytes(runner_rel,pins[runner_rel])
 import ast
 runner_ast=ast.parse(runner_data)
 lock_values=[ast.literal_eval(n.value) for n in runner_ast.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PARENT_LOCK' for t in n.targets)]
 if len(lock_values)!=1: raise SystemExit('portable pair lock constant unavailable')
 lock_pin=dict(path=str(pair_root/'inputs/launch_lock.json'),sha256=lock_values[0]); add_pin(lock_pin)
 lock_data,_=exact_bytes(rel(lock_pin['path']),lock_pin['sha256']); lock=json.loads(lock_data)
 # run_pair.read_contract checks these four inputs, every runtime file,
 # its external tax driver, and the complete source inventory.
 for name,key in (('inputs/objective.json','objective_sha256'),('inputs/proposal_bank.json','proposal_bank_sha256'),('inputs/source_manifest.json','source_manifest_sha256'),('ancestor_commute.py','ancestor_sha256')):
  add_pin(dict(path=str(pair_root/name),sha256=lock[key]))
 for name,digest in lock['runtime_file_sha256'].items(): add_pin(dict(path=str(pair_root/name),sha256=digest))
 add_pin(dict(path=str(pair_root.parent/'paygo_tax_comparison_20260924/run_paygo_two_rate.py'),sha256=lock['tax_driver_sha256']))
 pair_manifest_rel=rel(pair_root/'inputs/source_manifest.json')
 pair_manifest_data,_=exact_bytes(pair_manifest_rel,pins[pair_manifest_rel]); pair_manifest=json.loads(pair_manifest_data)
 pair_files=pair_manifest['files']
 if not (isinstance(pair_files,dict) and all(isinstance(k,str) and isinstance(v,str) and len(v)==64 for k,v in pair_files.items())): raise SystemExit('invalid pair source manifest: '+pair_manifest_rel)
 source_maps[pair_manifest_rel]=dict(source_root=rel(pair_root/'source'),files=pair_files)
 # Runtime source authentication validates only the explicit inventories named
 # by the controller and native ancestry contracts; do not invent a recursive
 # provenance closure from unrelated historical report metadata.
 for p in list(pins):
  if not p.endswith('.json'): continue
  data,_=exact_bytes(p,pins[p])
  doc=json.loads(data)
  pin=doc.get('source_manifest')
  if isinstance(pin,dict): add_pin(pin)
 for p in list(pins):
  if not p.endswith('.json'): continue
  data,_=exact_bytes(p,pins[p]); doc=json.loads(data)
  manifest_pin=doc.get('source_manifest'); source_root=doc.get('source_root')
  if not (isinstance(manifest_pin,dict) and isinstance(source_root,str)): continue
  manifest_rel=rel(manifest_pin['path']); manifest_data,_=exact_bytes(manifest_rel,pins[manifest_rel])
  files=json.loads(manifest_data)['files']
  if not (isinstance(files,dict) and all(isinstance(k,str) and isinstance(v,str) and len(v)==64 for k,v in files.items())): raise SystemExit('invalid source manifest: '+manifest_rel)
  root_rel=rel(source_root)
  old=source_maps.setdefault(manifest_rel,dict(source_root=root_rel,files=files))
  if old['source_root']!=root_rel or old['files']!=files: raise SystemExit('ambiguous source manifest placement: '+manifest_rel)
 return pins,source_maps
def native_reference_case_pins(manifest):
 """The three reference-case reads in native.setup, separate from export data."""
 contract_rel=rel(manifest['contract']['path'])
 contract_data,_=exact_bytes(contract_rel,manifest['contract']['sha256'])
 case=rel(json.loads(contract_data)['reference_case'])
 checkpoint_sha='090c9ebda662bf7837c4f4cf1d816159bc9d203a9babe7c70d00d0c4be1e575e'
 checkpoint=case+'/initial_state.pkl.gz'
 data,_=exact_bytes(checkpoint,checkpoint_sha)
 receipt=case+'/receipt.json'; parameters=case+'/parameters.csv'
 for path in (receipt,parameters):
  if not (ROOT/path).is_file(): raise SystemExit('native reference artifact absent: '+path)
 if json.loads((ROOT/receipt).read_text())['case_checkpoint_sha256']!=checkpoint_sha:
  raise SystemExit('native reference receipt/checkpoint mismatch')
 # setup pins the checkpoint itself. The receipt and parameter table are
 # required reads; pin their retained exact bytes before staging as well.
 return {checkpoint:checkpoint_sha,receipt:sha(ROOT/receipt),parameters:sha(ROOT/parameters)}
def preserve_attempt1():
 attempt1=DEPLOY/'attempt1'
 if attempt1.exists(): return
 attempt1.mkdir(parents=True)
 for name in ('inventory.json','stage.tar.gz','stage_receipt.json'):
  p=DEPLOY/name
  if p.is_file(): shutil.copy2(p,attempt1/name)
 resume=ROOT/'tmp/ces_share_overnight/stage_resume.log'
 if resume.is_file(): shutil.copy2(resume,attempt1/'stage_resume.log')
def main():
 if (OUT/'stage.tar.gz').exists(): raise SystemExit('Retained immutable attempt exists; select a fresh attempt before building')
 if sha(PARENT)!=PARENT_SHA: raise SystemExit('authenticated parent archive drift')
 with tarfile.open(PARENT) as a:
  prior=json.load(a.extractfile('inventory.json')); source={n.removeprefix('source/'):a.extractfile(n).read() for n in a.getnames() if n.startswith('source/')}
 if {k:sha(v) for k,v in source.items()}!=prior['files']: raise SystemExit('parent source-inventory drift')
 # Current executable production is preserved; frozen dependencies below replace only when pin-exact.
 for relative in ('code/model/production','code/model/production/reference_inputs','code/model/experiments/ces_normalized_shares','output/model/fixed_reference_economics_20260928/sources/fixed_price_v1'):
  add_tree(source,relative)
 for relative in (PLAN,ANCHOR,'output/model/publication_refactor_20260929/local_export_v1/inputs/arrays.npz','output/model/publication_refactor_20260929/local_export_v1/inputs/bundle.json',MANIFEST,OVERLAY+'/bootstrap.py',OVERLAY+'/frozen_sources/e5f_exact_policy_cache.py',OVERLAY+'/frozen_sources/test_e5f_exact_policy_cache.py'):
  p=ROOT/relative
  if not p.is_file(): raise SystemExit('required current/reference input absent: '+relative)
  source[relative]=p.read_bytes()
 manifest=json.loads(source[MANIFEST]); pins,source_maps=pin_dependencies(manifest)
 origins={}
 for relative,digest in sorted(pins.items()):
  data,origin=exact_bytes(relative,digest)
  source[relative]=data; origins[origin]=origins.get(origin,0)+1
 source_pins={}
 for record in source_maps.values():
  for relative,digest in record['files'].items():
   destination=str(Path(record['source_root'])/relative)
   old=source_pins.setdefault(destination,(relative,digest))
   if old!=(relative,digest): raise SystemExit('conflicting source-manifest destination: '+destination)
 for destination,(relative,digest) in sorted(source_pins.items()):
  data,origin=exact_bytes(destination,digest); source[destination]=data; origins[origin]=origins.get(origin,0)+1
 # The export is a declared artifact closure, including initial_state.pkl.gz separately.
 export=rel(manifest['local_export'])
 artifact_pins={f'{export}/{name}':digest for name,digest in manifest['artifact_hashes'].items()}
 artifact_pins[f'{export}/initial_state.pkl.gz']=manifest['checkpoint']['sha256']
 artifact_pins.update(native_reference_case_pins(manifest))
 for relative,digest in sorted(artifact_pins.items()):
  data,origin=exact_bytes(relative,digest); source[relative]=data; origins[origin]=origins.get(origin,0)+1
 for relative,digest in pins.items():
  if sha(source[relative])!=digest: raise SystemExit('staged declared pin mismatch: '+relative)
 for destination,(_,digest) in source_pins.items():
  if sha(source[destination])!=digest: raise SystemExit('staged declared source mismatch: '+destination)
 for relative,digest in artifact_pins.items():
  if sha(source[relative])!=digest: raise SystemExit('staged declared artifact mismatch: '+relative)
 plan=json.loads(source[PLAN]); anchor=json.loads(source[ANCHOR])
 if anchor.get('status')!='selected_numerically_verified' or anchor.get('chain')!=13 or anchor.get('arm')!='alternative': raise SystemExit('anchor completed.json authentication failed')
 if plan.get('selected_source_sha256')!=sha(source[ANCHOR]) or plan.get('source_checkpoint_sha256')!=sha(source[ANCHOR]): raise SystemExit('plan/anchor hash mismatch')
 rows=plan.get('target_contract',[])
 if len(rows)!=14 or sum(r.get('role')=='scored' for r in rows)!=11 or plan.get('target_fingerprint') is None or plan.get('weight_fingerprint') is None: raise SystemExit('target-row authentication failed')
 if len(plan.get('starts',()))!=4 or len(plan.get('coordinates',()))!=11 or any(plan.get('bounds',{}).get(k)!=[0.,.25] for k in ('delta_alpha_jump','delta_alpha')) or 'h_P' in plan.get('coordinates',()): raise SystemExit('CES plan drift')
 entry={p.name:sha(p) for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py','.sh') and p.name!='build_stage.py'}
 declared=set(pins)|set(source_pins)|set(artifact_pins); added=declared-set(prior['files']); added_bytes=sum(len(source[p]) for p in added)
 source_map_summary={k:dict(source_root=v['source_root'],file_count=len(v['files'])) for k,v in source_maps.items()}
 inventory=dict(files={k:sha(v) for k,v in sorted(source.items())},entrypoints=entry,parent_archive_sha256=PARENT_SHA,parent_inventory_sha256=sha(json.dumps(prior,sort_keys=True,indent=2).encode()),start_plan_sha256=sha(source[PLAN]),selected_source_sha256=sha(source[ANCHOR]),target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],remote_root=REMOTE,source_prefix='source/',no_cache_or_results=True,authenticated_dependencies=dict(fixed_reference_manifest=MANIFEST,declared_pin_count=len(pins),declared_source_manifests=source_map_summary,declared_source_file_count=len(source_pins),declared_artifact_count=len(artifact_pins),added_declared_files=len(added),added_declared_bytes=added_bytes,origins=origins,unresolved_pointers=[]))
 preserve_attempt1(); OUT.mkdir(parents=True,exist_ok=True); (OUT/'inventory.json').write_text(json.dumps(inventory,indent=2,sort_keys=True)+'\n')
 entries={'source/'+k:v for k,v in source.items()}; entries['inventory.json']=(OUT/'inventory.json').read_bytes(); entries.update({n:(HERE/n).read_bytes() for n in entry})
 stage=OUT/'stage.tar.gz'
 with stage.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
  with tarfile.open(fileobj=gz,mode='w') as a:
   for n,b in sorted(entries.items()):
    z=tarfile.TarInfo(n); z.size=len(b); z.mtime=0; z.mode=0o755 if n.endswith('.sh') else 0o644; a.addfile(z,io.BytesIO(b))
 staged=ROOT/'tmp/ces_share_overnight/staged_v3'
 if staged.exists(): shutil.rmtree(staged)
 staged.mkdir(parents=True)
 with tarfile.open(stage) as a: a.extractall(staged)
 receipt=dict(status='prepared_no_submission',archive=str(stage),sha256=sha(stage),bytes=stage.stat().st_size,source_files=len(source),staged_local_root=str(staged),target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],selected_source_sha256=sha(source[ANCHOR]),remote_root=REMOTE,added_declared_files=len(added),added_declared_bytes=added_bytes,unresolved_pointers=[])
 (OUT/'stage_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n'); print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__': main()
