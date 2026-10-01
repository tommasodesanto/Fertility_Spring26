"""Authenticate exact winner, model sources and credit arms without a lifecycle solve."""
import argparse,json,os,sys,time
from pathlib import Path
p=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(p))
import fixed_price_responses as driver
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
args.out.mkdir(parents=True,exist_ok=False)
a=driver.authenticate_candidate(args.out/'runtime')
P=a['P'];C=a['natural']
assert len(a['actual_parameters'])==len(a['params_rows'])==31
assert len(a['grid'])==120 and int(P.Nz)==9
assert P.unsecured_credit_limit==0. and P.native_due_stayer_credit and not bool(getattr(P,'native_solvency_credit',False))
assert C.unsecured_credit_limit is None and not C.native_due_stayer_credit and C.native_solvency_credit
assert P.hbar_first_child_jump==2.3 and P.psi_child==driver.BINDING['candidate_psi']
r=dict(status='passed_actual_runtime_zero_lifecycle',lifecycle_solves=0,actual_parameters=a['actual_parameters'],entry=a['entry'],driver_sha256=driver.sha(p/'fixed_price_responses.py'),source_binding_sha256=driver.sha(p/'source_binding.json'),candidate=driver.CONTRACT['candidate_case'],pid=os.getpid(),checked_epoch=time.time(),natural_support_certified=False)
driver.write(args.out/'initializer_receipt.json',r)
print(json.dumps(dict(status=r['status'],candidate=r['candidate'],lifecycle_solves=0)))
