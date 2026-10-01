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
 if 'prepared_native_inputs' in plan:extra.append(Path(plan['prepared_native_inputs']['path']))
 if 'prepared_consumer_compatibility' in plan:
  approved=Path(plan['prepared_consumer_compatibility']['path']);extra.append(approved)
  receipt=json.loads(approved.read_text())
  extra += [Path(receipt[k]['path']) for k in ('original_generator_snapshot','reviewed_audit','reviewed_diff')]
 if 'diagnostic_measurement_reuse' in plan:
  manifest=Path(plan['diagnostic_measurement_reuse']['path']);extra.append(manifest)
  def referenced_pins(value):
   if isinstance(value,dict):
    if set(value)=={'path','sha256'}:
     pin=Path(value['path'])
     if pin.is_relative_to(ROOT):extra.append(pin)
    else:
     for item in value.values():referenced_pins(item)
   elif isinstance(value,list):
    for item in value:referenced_pins(item)
  referenced_pins(json.loads(manifest.read_text()))
extra.append(OUT/'fit_preflight.py')
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
