"""Torch-only inspection of the failed saved point; no solves or mutations."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
import sys,gzip,pickle
from pathlib import Path
import numpy as np
import driver
BASE=Path(__file__).resolve().parent
OUT=BASE/'postrun_audit_v1'
OUT.mkdir(exist_ok=False)
smoke=driver.read(BASE/'smoke_v3/worker/receipt.json')
for name,digest in smoke['source_pins'].items(): assert driver.sha(BASE/name)==digest
c,obj,evaluator,reference,point,ref_receipt=driver.load_runtime(OUT)
manifest=driver.install_overlay(OUT,evaluator)
case=BASE/'run_v1/worker/case'
receipt=driver.read(case/'receipt.json')
assert driver.sha(case/'initial_state.pkl.gz')==receipt['case_checkpoint_sha256']
with gzip.open(case/'initial_state.pkl.gz','rb') as f: packet=pickle.load(f)
P=packet['parameters'];e=packet['evaluation'];sol=packet['solution'];model=evaluator.rt['model']
arrays=[np.asarray(sol.fert_extra_probs),np.asarray(e.policy.fert_extra_probs),np.asarray(P._fert_extra_probs)]
expected=np.shape(e.policy.fert_probs)[:-1]+(2,2,int(P.n_child_states))
checks=dict(shapes_match=all(a.shape==expected for a in arrays),finite=all(np.isfinite(a).all() for a in arrays),in_unit_interval=all(a.min()>=0 and a.max()<=1 for a in arrays),equal=all(np.array_equal(a,arrays[0]) for a in arrays),independent=all(not np.shares_memory(arrays[i],arrays[k]) for i in range(3) for k in range(i)))
sums=arrays[0].sum(axis=5);fec=model.get_fecundity_by_age(P);rows=[]
for j in range(int(P.A_f_start)-1,int(P.A_f_end)):
 for n in range(2):
  for m in range(n+1):
   total=sums[:,:,:,j,:,n,m]
   err=np.minimum(abs(total),abs(total-1.))
   outer=e.policy.fert_probs[:,:,:,j,:,1] if n==0 else e.policy.fert2_probs[:,:,:,j,:,1,n-1,m]
   reached=e.g_pre[:,:,:,j,:,n,m]*outer*fec[j]
   zero=total==0.;bad=(abs(total-1.)>1e-12)&~zero
   idx=np.unravel_index(np.argmax(err),err.shape);b,t,l,z=map(int,idx)
   reached_error=np.where(reached>1e-14,err,0.)
   ridx=np.unravel_index(np.argmax(reached_error),err.shape);rb,rt,rl,rz=map(int,ridx)
   rows.append(dict(age=float(P.age_start+P.da*j),j=j,n=n,m=m,maximum_sum_error=float(err.max()),maximum_reached_sum_error=float(reached_error.max()),reached_total=float(reached.sum()),reached_zero_menu_mass=float(reached[zero].sum()),reached_nonzero_bad_sum_mass=float(reached[bad].sum()),nonzero_bad_count=int(bad.sum()),argmax=[b,t,l,z],argmax_probabilities=arrays[0][b,t,l,j,z,:,n,m].tolist(),argmax_reached=float(reached[idx]),reached_argmax=[rb,rt,rl,rz],reached_argmax_probabilities=arrays[0][rb,rt,rl,j,rz,:,n,m].tolist(),reached_argmax_mass=float(reached[ridx])))
checks.update(maximum_sum_error=max(r['maximum_sum_error'] for r in rows),reached_zero_menu_mass=sum(r['reached_zero_menu_mass'] for r in rows),reached_nonzero_bad_sum_mass=sum(r['reached_nonzero_bad_sum_mass'] for r in rows))
result=dict(status='inspection_only_original_gate_failed',model_solves=0,checks=checks,rows=rows,case_checkpoint_sha256=receipt['case_checkpoint_sha256'],original_receipt_sha256=driver.sha(case/'receipt.json'),source_pins=smoke['source_pins'],effective_source_manifest_sha256=driver.sha(OUT/'effective_source_manifest.json'),original_completion=driver.read(BASE/'run_v1/completion.json'))
driver.write(OUT/'audit.json',result)
driver.write(OUT/'case_artifact_hashes.json',{str(p.relative_to(case)):driver.sha(p) for p in sorted(case.rglob('*')) if p.is_file()})
print(checks,flush=True)
