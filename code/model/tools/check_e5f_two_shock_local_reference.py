"""Recompute pristine parent baseline in the same runtime as the new experiment."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ.setdefault(key,'1')
import argparse,sys,json,time
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--checkpoint',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
sys.path[:0]=[str(args.source/'code/model'),str(args.source/'code/model/tools')]
import numpy as np
import run_e5f_independent_numerical_audit as audit
from intergen_eqscale_seq_optimized import solver as model
packet=audit.load_checkpoint(args.checkpoint);P=packet['parameters'];bg=packet['b_grid'];ref=packet['evaluation'].policy;price=ref.price
P.joint_nested_choice=False
expected_second=P._fert2_probs.copy();start=time.monotonic()
result=model.solve_bellman_full_markov_income(P.user_cost_rate*price,price,P,bg,model.precompute_shared(P,bg),continuation_V=ref.V)
fields=('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value')
arrays=dict(zip(fields,result[:9]));arrays['fert2_probs']=P._fert2_probs
comparison={key:dict(exact=bool(np.array_equal(value,expected_second if key=='fert2_probs' else getattr(ref,key))),max_abs=float(np.max(abs(value-(expected_second if key=='fert2_probs' else getattr(ref,key)))))) for key,value in arrays.items()}
args.output.mkdir(parents=True,exist_ok=True)
np.savez_compressed(args.output/'reference_arrays.npz',**arrays)
audit.save_json(args.output/'reference.json',dict(source=str(args.source),elapsed_seconds=time.monotonic()-start,comparison=comparison,source_hashes={str(p.relative_to(args.source)):audit.digest(p) for p in (args.source/'code/model/intergen_eqscale_seq_optimized').glob('*.py')}))
print(json.dumps(comparison),flush=True)
