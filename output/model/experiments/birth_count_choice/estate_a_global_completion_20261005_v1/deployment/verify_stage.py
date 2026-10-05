"""Check the inherited frozen model overlay and this stage's entrypoint pins."""
import argparse,hashlib,json
from pathlib import Path
REMOTE=Path('/scratch/td2248/projects/estate_birth_global_completion_20261005_v1')
REPO=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--container',action='store_true');p.add_argument('--host',action='store_true');a=p.parse_args()
 assert a.container!=a.host
 root=Path('/work/deployment') if a.container else REMOTE
 m=json.loads((root/'stage_manifest.json').read_text())
 inv=json.loads((root/'inventory.json').read_text())
 assert sha(root/'control/plan.json')==m['completion_plan_sha256']
 assert sha(root/'control/original_plan.json')==m['original_plan_sha256']
 assert sha(root/'inventory.json')==m['parent_inventory_sha256']
 source=REPO if a.container else root/'source'
 for rel,digest in inv['files'].items():
  assert sha(source/rel)==digest,rel
 for rel,digest in m['entrypoints'].items():
  assert sha(root/rel)==digest,rel
 print(json.dumps(dict(status='passed_zero_solves',mode='container' if a.container else 'host',source_files=len(inv['files']),
                       target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'])))
if __name__=='__main__':main()
