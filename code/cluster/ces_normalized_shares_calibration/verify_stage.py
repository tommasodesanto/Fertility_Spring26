"""Verify every staged source byte and target contract; zero model solves."""
import hashlib,json,sys
from pathlib import Path
REMOTE=Path("/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3"); STAGE=Path("/work/deployment"); REPO=Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 mode=sys.argv[1:] or ["--host"]
 if mode not in (["--host"],["--container"],["--local"]):raise SystemExit("Usage: --host|--container|--local")
 if mode==["--local"]:
  stage=Path(__file__).resolve().parent; root=stage/"source"
 else:
  stage=STAGE if mode==["--container"] else REMOTE; root=REPO if mode==["--container"] else stage/"source"
 inv=json.loads((stage/"inventory.json").read_text())
 for rel,d in inv["files"].items():
  if not (root/rel).is_file() or sha(root/rel)!=d:raise SystemExit("Source pin drift: "+rel)
 for name,d in inv["entrypoints"].items():
  if sha(stage/name)!=d:raise SystemExit("Entrypoint pin drift: "+name)
 deps=inv.get("authenticated_dependencies",{})
 if deps.get("unresolved_pointers")!=[]:raise SystemExit("Unresolved authenticated pointers")
 for rel in (deps.get("fixed_reference_manifest"),):
  if not rel or rel not in inv["files"]:raise SystemExit("Fixed-reference manifest omitted")
 plan=json.loads((root/"output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json").read_text())
 if plan["target_fingerprint"]!=inv["target_fingerprint"] or plan["weight_fingerprint"]!=inv["weight_fingerprint"] or len(plan["starts"])!=4:raise SystemExit("Plan fingerprint/count mismatch")
 print(json.dumps(dict(status="passed",source_files=len(inv["files"]),target_fingerprint=inv["target_fingerprint"],weight_fingerprint=inv["weight_fingerprint"])))
if __name__=="__main__":main()
