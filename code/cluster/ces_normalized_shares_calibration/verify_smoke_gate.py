"""Fail-closed gate for the single chain-0 CES native smoke."""
import hashlib,json,sys
from pathlib import Path
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(root):
 root=Path(root); inv=read(root/"inventory.json"); plan=read(root/"source/output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json"); launch=root/"results/smoke_chain_0"; run=launch/"run"
 for p in (launch/"launcher_start.json",run/"start_contract.json",run/"cases.json",run/"heartbeat.json",run/"completed.json",run/"native_postcheck/completed.json"):
  if not p.is_file():raise AssertionError("Missing smoke receipt: "+str(p))
 start,contract,cases,heart,done,child=map(read,(launch/"launcher_start.json",run/"start_contract.json",run/"cases.json",run/"heartbeat.json",run/"completed.json",run/"native_postcheck/completed.json"))
 if start["wall_seconds"]!=5400 or start["maximum_objective_calls"]!=500 or start["final_native_reserve_seconds"]!=1800:raise AssertionError("smoke budget drift")
 if start["stage_inventory_sha256"]!=sha(root/"inventory.json") or len(cases)!=2 or heart["status"]!="completed":raise AssertionError("incomplete smoke")
 if done.get("status")!="selected_numerically_verified" or done.get("objective_calls")!=2:raise AssertionError("smoke did not complete exactly two verified calls")
 if len({tuple(sorted(case.get("parameters",{}).items())) for case in cases})!=2 or any(case.get("status")!="passed" for case in cases):raise AssertionError("smoke requires two distinct passed physical parameter cases")
 keys=("target_fingerprint","weight_fingerprint","starts_file_sha256","selected_source_sha256","source_checkpoint_sha256")
 for key in keys:
  expected=inv["start_plan_sha256"] if key=="starts_file_sha256" else plan.get(key)
  if contract.get(key)!=expected:raise AssertionError("start contract/plan mismatch: "+key)
 if contract["target_fingerprint"]!=inv["target_fingerprint"] or contract["weight_fingerprint"]!=inv["weight_fingerprint"] or contract["starts_file_sha256"]!=inv["start_plan_sha256"]:raise AssertionError("smoke inventory contract drift")
 exact=child.get("exact_tables",{})
 if child.get("status")!="full_native_postcheck_passed" or child.get("target_rows")!=14 or child.get("parameter_rows")!=31 or child.get("residual_count")!=11 or len(child.get("standard_plot_hashes",{}))!=17 or not exact.get("target_fit_csv_equal") or not exact.get("target_fit_experimental_csv_equal") or not exact.get("parameters_csv_equal"):raise AssertionError("fresh selected native postcheck shape/equality failed")
 for key in keys:
  expected=inv["start_plan_sha256"] if key=="starts_file_sha256" else plan.get(key)
  if child.get(key)!=expected:raise AssertionError("child/plan mismatch: "+key)
 return dict(status="smoke_gate_passed",chain=0,target_fingerprint=inv["target_fingerprint"],weight_fingerprint=inv["weight_fingerprint"],starts_file_sha256=inv["start_plan_sha256"])
if __name__=="__main__":print(json.dumps(verify(sys.argv[1] if len(sys.argv)==2 else "/scratch/td2248/projects/ces_normalized_shares_overnight_20261003_v3"),sort_keys=True))
