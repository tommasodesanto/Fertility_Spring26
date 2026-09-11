"""Hash original four-case outputs and package light artifacts; never solve."""
import hashlib,json,pathlib,tarfile,datetime
root=pathlib.Path('/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/output/utility_fiscal_decomposition_20260911a')
cases=('old_old','old_balanced','new_old','new_balanced')
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8*1024*1024),b''): h.update(chunk)
 return h.hexdigest()
manifest={'case_order':cases,'created_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':{},'checkpoints':{},'original_graphs':{}}
for name in cases:
 case=root/name
 summary=json.loads((case/'summary.json').read_text())
 assert summary['status']=='passed_initial_candidate_loop',name
 assert summary['repetitions']==summary['stationary_solves']==1 and summary['normalized'] is False,name
 assert summary['calibrated_smm'] is False and summary['perfect_foresight_solved'] is False,name
 checkpoint=case/'repetition_01/initial_state.pkl.gz'
 pin=digest(checkpoint)
 assert pin==summary['final']['checkpoint_sha256'],name
 manifest['checkpoints'][name]={'path':str(checkpoint),'sha256':pin,'bytes':checkpoint.stat().st_size,'downloaded':False}
 graphs=sorted((case/'repetition_01/standard_diagnostics').glob('*.png'))
 assert len(graphs)==17,(name,len(graphs))
 manifest['original_graphs'][name]=[{'path':str(p),'sha256':digest(p),'bytes':p.stat().st_size} for p in graphs]
 for p in sorted(case.rglob('*')):
  if p.is_file() and p.suffix in ('.json','.csv','.png'):
   manifest['files'][str(p.relative_to(root))]={'sha256':digest(p),'bytes':p.stat().st_size}
for p in sorted(root.glob('slurm_*')):
 if p.is_file(): manifest['files'][p.name]={'sha256':digest(p),'bytes':p.stat().st_size}
receipt=root/'collection_manifest.json'
assert not receipt.exists(),'collection already recorded'
receipt.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
archive=root/'light_outputs.tar.gz'
assert not archive.exists(),'archive already exists'
with tarfile.open(archive,'w:gz') as bundle:
 for rel in manifest['files']: bundle.add(root/rel,arcname=rel)
 bundle.add(receipt,arcname=receipt.name)
print(json.dumps({'archive_path':str(archive),'archive_sha256':digest(archive),'archive_bytes':archive.stat().st_size,'light_artifacts':len(manifest['files']),'original_graphs':68,'checkpoint_count':4}))
