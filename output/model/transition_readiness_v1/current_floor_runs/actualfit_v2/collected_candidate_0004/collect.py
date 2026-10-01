import hashlib,json,pathlib,subprocess
HERE=pathlib.Path(__file__).resolve().parent
REMOTE='/scratch/td2248/projects/transition_readiness_v1/current_floor/results/actualfit_v2'
code='''import pathlib,json,hashlib
r=pathlib.Path("REMOTE");c=r/"candidate_0004"
files=list(c.glob("*.json"))+list((c/"state_2023_checkpoint").glob("*.json"))
for h in (12,16):
 d=c/("horizon_%03d"%h);files+=list(d.glob("*.json"))
 for name in ["endpoint"]+[p.name for p in d.iterdir() if p.name.startswith("map_")]:
  files+=list((d/name).glob("*.json"))
 files+=list((d/"endpoint"/"one_step").glob("*.json"))
p=c/"state_2023_checkpoint"/"actual_2023.pkl.gz"
if p.stat().st_size<=250*1024**2:files.append(p)
print(json.dumps({str(p.relative_to(r)):{"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"bytes":p.stat().st_size} for p in files}))
'''.replace('REMOTE',REMOTE)
reply=subprocess.run(['ssh','-4','-o','BatchMode=yes','torch','python -c '+__import__('shlex').quote(code)],check=True,capture_output=True,text=True)
manifest=json.loads(reply.stdout)
(HERE/'remote_inventory.json').write_text(json.dumps(manifest,indent=2)+'\n')
for rel,pin in manifest.items():
 target=HERE/rel;target.parent.mkdir(parents=True,exist_ok=True)
 subprocess.run(['scp','-4','-o','BatchMode=yes','-q','torch:'+REMOTE+'/'+rel,str(target)],check=True)
 if hashlib.sha256(target.read_bytes()).hexdigest()!=pin['sha256']:raise RuntimeError('Copied bytes differ: '+rel)
receipt=HERE/'candidate_0004/state_2023_checkpoint/checkpoint_receipt.json'
assert hashlib.sha256(receipt.read_bytes()).hexdigest()=='a3e996e74c0de2d53b1e66e7f9de4901cb2cb6950a99ed2a9d723a9539db1de1'
print(json.dumps({'verified_copied_files':len(manifest),'bytes':sum(p['bytes'] for p in manifest.values()),'checkpoint_receipt_authoritative_hash_pass':True}))
