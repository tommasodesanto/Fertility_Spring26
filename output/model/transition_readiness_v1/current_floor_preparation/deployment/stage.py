#!/usr/bin/env python3
"""Deterministic own-source snapshot; never writes the calibration deployment."""
import hashlib,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
OUT=Path(__file__).resolve().parent
DEST=OUT/'stage'
REMOTE='/scratch/td2248/projects/transition_readiness_v1/current_floor'
code=ROOT/'code/model/experiments/transition_readiness'
handoff=ROOT/'output/model/transition_readiness_v1/current_floor_handoff'
h=json.loads((handoff/'handoff.json').read_text())
winner=ROOT/h['selected_point']['local_packet']
pins=ROOT/h['checkpoint_and_sources']['source_pins_manifest']
paths=[p for p in code.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
paths += [p for p in handoff.rglob('*') if p.is_file()]
paths += [p for p in winner.rglob('*') if p.is_file()]
paths += [pins]
planpaths=[p for p in OUT.glob('*plan.json') if p.is_file()]
paths += planpaths
extra=[]
for planpath in planpaths:
 plan=json.loads(planpath.read_text())
 extra += [Path(plan['target_contract'][k]['path']) for k in ('blocks','annual')]
paths += extra
if DEST.exists():shutil.rmtree(DEST)
files={}
for p in sorted(set(paths)):
 rel=str(p.relative_to(ROOT));target=DEST/'source'/rel
 target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 files[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
(DEST/'inventory.json').write_text(json.dumps(dict(schema='current_floor_own_overlay_v1',files=files),indent=2)+'\n')
# Bind compact directories so overlaying nonexistent new file paths works reliably.
mounts=[str(p.relative_to(ROOT)) for p in (code,handoff,winner)]
mounts += [str(pins.relative_to(ROOT))]
mounts += [str(p.relative_to(ROOT)) for p in planpaths+extra]
(DEST/'mounts.txt').write_text('\n'.join(mounts)+'\n')
shutil.copy2(code/'floor_launch.sh',DEST/'floor_launch.sh')
(OUT/'stage_receipt.json').write_text(json.dumps(dict(remote=REMOTE,file_count=len(files),mounts=mounts,inventory_sha256=hashlib.sha256((DEST/'inventory.json').read_bytes()).hexdigest(),files=files),indent=2)+'\n')
print(json.dumps(dict(stage=str(DEST),files=len(files),mounts=len(mounts))))
