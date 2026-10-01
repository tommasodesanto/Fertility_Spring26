"""Compact frozen-source-only deployment; never include arrays or solution caches."""
from pathlib import Path
import gzip,hashlib,io,json,tarfile
HERE=Path(__file__).resolve().parent;PACKET=HERE.parents[1];ROOT=PACKET.parents[4]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
binding=json.loads((PACKET/'source_binding.json').read_text());old=json.loads((ROOT/binding['source_pins_path']).read_text());pins={r['path']:r['sha256'] for r in binding['files']}|old
# Complete the first-import package using only byte-identical original pinned indexed sources.
origin='output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source/small_credit_lab/'
target='output/model/fixed_reference_economics_20260928/credit_no_taper_v1/small_credit_v1/source/small_credit_lab/'
package_receipt=[]
for rel,digest in sorted(old.items()):
 if rel.startswith(origin) and rel.endswith('.py'):
  mapped=target+rel[len(origin):]
  assert sha(ROOT/rel)==sha(ROOT/mapped)==digest,mapped
  pins[mapped]=digest
  package_receipt.append(dict(source=rel,target=mapped,sha256=digest))
(HERE/'package_completeness.json').write_text(json.dumps(package_receipt,indent=2)+'\n')
sources={};large={};base='/scratch/td2248/projects/grid_resolution_credit053_v2'
for rel,digest in pins.items():
 p=ROOT/rel;assert sha(p)==digest,rel
 if p.suffix=='.npz':large[rel]=dict(sha256=digest,remote_path=base+'/source/'+rel);continue
 assert p.suffix in ('.py','.json','.csv','.sh','.md') and p.stat().st_size<2_000_000,rel
 sources[rel]=digest
for p in [PACKET/'source_binding.json',ROOT/binding['source_pins_path'],PACKET/'deployment/initialize_runtime.py']+list(HERE.glob('*.py'))+list(HERE.glob('*.sh')):sources[str(p.relative_to(ROOT))]=sha(p)
inventory=dict(files=sources,reused_large_inputs=large,source_prefix='source/',includes_arrays=False,includes_solution_caches=False,base_remote=base,frozen_driver_sha256=sha(PACKET/'fixed_price_responses.py'),binding_sha256=sha(PACKET/'source_binding.json'))
assert inventory['frozen_driver_sha256']=='7bed458277ff0bfb3923301f1f8b8907ac01cb5ebcbc7002cab4e5b47a432f0f'
assert inventory['binding_sha256']=='b215b7332d03e163244808cd294c17ee1d50b513139dd5c83385f3eb95ec04fe'
(HERE/'inventory.json').write_text(json.dumps(inventory,indent=2,sort_keys=True)+'\n')
# The utility-floor psi packet directory is bound once; other files bind individually.
prefix='output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/'
(HERE/'mounts.txt').write_text('\n'.join(rel for rel in sorted(sources) if not rel.startswith(prefix))+'\n')
archive=HERE/'mechanism_stage.tar.gz'
with archive.open('wb') as raw:
 with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as zipped:
  with tarfile.open(fileobj=zipped,mode='w') as tar:
   entries=[('source/'+rel,ROOT/rel) for rel in sorted(sources)]+[(name,HERE/name) for name in ('inventory.json','mounts.txt','launch_torch.sh','submit_torch.sh')]
   for name,p in entries:
    b=p.read_bytes();i=tarfile.TarInfo(name);i.size=len(b);i.mtime=0;i.mode=0o755 if p.suffix=='.sh' else 0o644;tar.addfile(i,io.BytesIO(b))
r=dict(archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,compact_sources=len(sources),reused_inputs=large,includes_arrays=False,includes_solution_caches=False)
(HERE/'stage_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
