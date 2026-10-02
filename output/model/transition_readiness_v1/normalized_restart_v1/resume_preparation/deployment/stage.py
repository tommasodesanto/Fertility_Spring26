import hashlib,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[6];OUT=Path(__file__).resolve().parent;STAGE=OUT/'stage';code=ROOT/'code/model/experiments/transition_readiness';packet=OUT.parents[1]
paths=[p for p in code.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
paths += [packet/'handoff.json',packet/'source_pins.json']
paths += [p for folder in ('native_reports','evidence') for p in (packet/folder).rglob('*') if p.is_file()]
extra=[OUT/'plan_preflight.py']+list(OUT.glob('*plan.json'))
for planpath in OUT.glob('*plan.json'):
 plan=json.loads(planpath.read_text());extra += [Path(plan['target_contract'][k]['path']) for k in ('blocks','annual')]
 if 'prepared_native_inputs' in plan:extra.append(Path(plan['prepared_native_inputs']['path']))
 if 'prepared_consumer_compatibility' in plan:
  p=Path(plan['prepared_consumer_compatibility']['path']);extra.append(p);r=json.loads(p.read_text())
  extra += [Path(r[k]['path']) for k in ('original_generator_snapshot','original_execution_plan','reviewed_audit','reviewed_diff')]
paths += extra
if STAGE.exists():shutil.rmtree(STAGE)
files={}
for p in sorted(set(paths)):
 rel=str(p.relative_to(ROOT));t=STAGE/'source'/rel;t.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(OUT/'floor_launch.sh' if p==code/'floor_launch.sh' else p,t);files[rel]=hashlib.sha256(t.read_bytes()).hexdigest()
(STAGE/'inventory.json').write_text(json.dumps({'schema':'normalized_resume_own_overlay_v1','files':files},indent=2)+'\n')
mounts=[str(code.relative_to(ROOT)),str(packet.relative_to(ROOT))]+[str(p.relative_to(ROOT)) for p in extra if not p.is_relative_to(packet)]
(STAGE/'mounts.txt').write_text('\n'.join(mounts)+'\n');shutil.copy2(OUT/'floor_launch.sh',STAGE/'floor_launch.sh')
print(json.dumps({'files':len(files),'inventory_sha256':hashlib.sha256((STAGE/'inventory.json').read_bytes()).hexdigest()}))
