"""Zero-lifecycle GE source, winner and fixed-q0 receipt initializer."""
import argparse,json,sys,time
from pathlib import Path
p=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(p))
import run_ge as driver
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--q0-credit-case',type=Path,required=True);a=ap.parse_args()
a.out.mkdir(parents=True,exist_ok=False)
plan=driver.verify_plan(p/'plan.json',zeroLC=True)
_,auth=driver.candidate_runtime(a.out/'runtime')
assert len(auth['actual_parameters'])==31 and len(auth['grid'])==120 and auth['natural'].Nz==9
assert (a.q0_credit_case/'receipt.json').is_file() and (a.q0_credit_case/'closure.json').is_file()
r=dict(status='passed_zero_lifecycle_initializer',lifecycle_solves=0,candidate=plan['candidate_case'],plan_sha256=driver.sha(p/'plan.json'),q0_receipt_sha256=driver.sha(a.q0_credit_case/'receipt.json'),checked_epoch=time.time())
driver.write(a.out/'initializer_receipt.json',r)
print(json.dumps(r))
