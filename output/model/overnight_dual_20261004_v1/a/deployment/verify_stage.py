import argparse,hashlib,json
from pathlib import Path
ROOT=Path('/scratch/td2248/projects/estate_birth_overnight_20261004_v1');REPO=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--host",action="store_true");p.add_argument("--container",action="store_true");a=p.parse_args();assert a.host!=a.container
 root=Path("/work/deployment") if a.container else ROOT;m=json.loads((root/"manifest.json").read_text());inv=json.loads((root/"inventory.json").read_text());source=REPO if a.container else root/"source"
 assert sha(root/"inventory.json")==m["inventory_sha256"]
 for rel,digest in inv["files"].items():assert sha(source/rel)==digest,rel
 for rel,digest in m["entrypoints"].items():assert sha(root/rel)==digest,rel
 print(json.dumps(dict(status="passed_zero_solves",arm='a',files=len(inv["files"]),target=inv["target_fingerprint"],weight=inv["weight_fingerprint"])))
