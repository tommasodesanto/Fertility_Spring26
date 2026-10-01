"""Freeze complete previously verified dependencies and this isolated driver."""
from pathlib import Path
import gzip,hashlib,io,json,tarfile
HERE=Path(__file__).resolve().parent;PACKET=HERE.parent;ROOT=PACKET.parents[3]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
ECON=ROOT/'output/model/fixed_reference_economics_20260928'
old=json.loads((ECON/'utility_floor_psi_v1/mechanism_responses_v1/deployment/replacement_v2/inventory.json').read_text())
round3=json.loads((ECON/'utility_floor_psi_round3_v1/deployment/inventory.json').read_text())
pins=dict(old['files'])
for rel,digest in round3['files'].items():
 if rel in pins:assert pins[rel]==digest,('conflicting reviewed dependencies',rel)
 pins[rel]=digest
for folder in ('entry_calibration_pilot_v1','parenthood_floor_quick_v1'):
 for p in (ECON/folder).iterdir():
  if p.suffix in ('.py','.json') and p.is_file():pins[str(p.relative_to(ROOT))]=sha(p)
ref=ECON/'utility_calibration_round1_v1/deployment/attempt2/comparison_requested'
for name in ('with_A_parameters.csv','with_A_target_fit.csv'):pins[str((ref/name).relative_to(ROOT))]=sha(ref/name)
for p in list(PACKET.glob('*.py'))+list(PACKET.glob('*.json'))+list(PACKET.glob('*.md'))+list((PACKET/'reference').glob('*'))+list(HERE.glob('*.py'))+list(HERE.glob('*.sh')):
 pins[str(p.relative_to(ROOT))]=sha(p)
for rel,digest in pins.items():assert sha(ROOT/rel)==digest,('source drift',rel)
manifest=dict(files=pins,reused_large_inputs=old['reused_large_inputs'],source_prefix='source/',includes_arrays=False,includes_solution_caches=False,dependency_provenance='round3 plus reviewed repaired prescribed-price full small_credit_lab package')
(HERE/'inventory.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
econprefix='output/model/fixed_reference_economics_20260928/utility_share_A_decomposition_v1/'
(HERE/'mounts.txt').write_text('\n'.join(rel for rel in sorted(pins) if not rel.startswith(econprefix))+'\n')
archive=HERE/'stage.tar.gz'
with archive.open('wb') as raw:
 with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as zipped:
  with tarfile.open(fileobj=zipped,mode='w') as tar:
   for name,p in [('source/'+rel,ROOT/rel) for rel in sorted(pins)]+[(name,HERE/name) for name in ('inventory.json','mounts.txt','launch_torch.sh','submit_torch.sh')]:
    data=p.read_bytes();i=tarfile.TarInfo(name);i.size=len(data);i.mtime=0;i.mode=0o755 if p.suffix=='.sh' else 0o644;tar.addfile(i,io.BytesIO(data))
r=dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,files=len(pins),reused_large_inputs=old['reused_large_inputs'])
(HERE/'stage_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
