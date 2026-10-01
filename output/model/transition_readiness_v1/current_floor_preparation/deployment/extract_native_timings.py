import gzip,hashlib,json,pickle,sys
from pathlib import Path
folder=Path('/scratch/td2248/projects/transition_readiness_v1/current_floor/results/smoke_v3/native_reference/repeat_0')
receipt=folder/'native_solve_unverified.json';expected_receipt='ab223ec92324949c795bb8f38f3680a82d59a2ed3a465733a130a25682ec95e6'
assert hashlib.sha256(receipt.read_bytes()).hexdigest()==expected_receipt
r=json.loads(receipt.read_text());checkpoint=folder/'native_solve_unverified.pkl.gz'
assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()==r['checkpoint']['sha256']=='5dcdb6a3edb31e299124c65f671fc876c1693e0f997e2ab2c149dd23bfef3d6b'
with gzip.open(checkpoint,'rb') as f:packet=pickle.load(f)
solution=packet['solution'];timings=getattr(solution,'timings',None)
out=dict(status='authenticated_saved_solution_timings',source_receipt_path=str(receipt),source_receipt_sha256=expected_receipt,source_checkpoint_path=str(checkpoint),source_checkpoint_sha256=r['checkpoint']['sha256'],original_native_policy_calls=r['policy_calls'],new_native_policy_calls=0,solution_type=type(solution).__name__,timings=timings,additional_timing_fields={k:v for k,v in vars(solution).items() if 'timing' in k.lower() or k in ('n_full','n_fast')})
def plain(x):
 if hasattr(x,'tolist'):return plain(x.tolist())
 if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [plain(v) for v in x]
 return x
print(json.dumps(plain(out),indent=2))
