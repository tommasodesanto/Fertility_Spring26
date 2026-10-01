import hashlib,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5];OUT=Path(__file__).resolve().parent;STAGE=OUT/'stage';code=ROOT/'code/model/experiments/transition_readiness';packet=OUT.parent
paths=[p for p in code.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
paths += [p for p in packet.rglob('*') if p.is_file() and 'deployment' not in p.relative_to(packet).parts]
paths += [OUT/'plan_preflight.py']+list(OUT.glob('*plan.json'))
for p in OUT.glob('*plan.json'):
 plan=json.loads(p.read_text());paths += [Path(plan['target_contract'][k]['path']) for k in ('blocks','annual')]
if STAGE.exists():shutil.rmtree(STAGE)
files={}
for p in sorted(set(paths)):
 rel=str(p.relative_to(ROOT));target=STAGE/'source'/rel;target.parent.mkdir(parents=True,exist_ok=True)
 shutil.copy2(OUT/'floor_launch.sh' if p==code/'floor_launch.sh' else p,target);files[rel]=hashlib.sha256(target.read_bytes()).hexdigest()
(STAGE/'inventory.json').write_text(json.dumps({'schema':'normalized_restart_own_overlay_v1','files':files},indent=2)+'\n')
mounts=[str(code.relative_to(ROOT)),str(packet.relative_to(ROOT)),str((OUT/'plan_preflight.py').relative_to(ROOT))]
# The whole packet bind includes staged plan and its compact evidence.
for p in OUT.glob('*plan.json'):
 plan=json.loads(p.read_text());mounts += [str(Path(plan['target_contract'][k]['path']).relative_to(ROOT)) for k in ('blocks','annual')]
(STAGE/'mounts.txt').write_text('\n'.join(mounts)+'\n');shutil.copy2(OUT/'floor_launch.sh',STAGE/'floor_launch.sh')
print(json.dumps({'files':len(files),'inventory_sha256':hashlib.sha256((STAGE/'inventory.json').read_bytes()).hexdigest()}))
