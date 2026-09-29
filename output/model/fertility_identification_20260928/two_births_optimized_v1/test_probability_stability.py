"""Torch-only proposed numerical repair test; no model imports or solves.

Preserve the failed experiment. This creates a distinct transformed source
candidate and tests its probability calculation; it does not install or adopt it.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
import sys,ast,json,hashlib,unittest
from pathlib import Path
assert sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit()
import numpy as np
import patch_solver
import tests_solver
BASE=Path(__file__).resolve().parent
OUT=BASE/'probability_repair_v1';OUT.mkdir(exist_ok=False)
old='    inclusive, probabilities = logsumexp(menu / kappa, axis=menu.ndim - 1)'
new="""    scaled = menu / kappa
    inclusive, _ = logsumexp(scaled, axis=menu.ndim - 1)
    shifted_weights = np.exp(scaled - np.max(scaled, axis=-1, keepdims=True))
    probabilities = shifted_weights / shifted_weights.sum(axis=-1, keepdims=True)"""
original_patch=patch_solver.patch_text
failed=original_patch(tests_solver.ORIGINAL)
assert failed.count(old)==1
corrected=failed.replace(old,new,1)
compile(corrected,'<proposed stable source>','exec')
(OUT/'solver_candidate.py').write_text(corrected)
namespace=dict(tests_solver.NS)
exec(compile(tests_solver.functions(corrected,['two_birth_last_opportunity']),'<stable helper>','exec'),namespace)
oldfn=tests_solver.NS['two_birth_last_opportunity'];newfn=namespace['two_birth_last_opportunity']
wait=np.array([-8e8,-1e8,-1e6,-10.,-2e9]);success=wait+np.array([0.,5.,-7.,2.,0.])
values0,p0=oldfn(wait,success,.83,.332);values1,p1=newfn(wait,success,.83,.332)
live=wait>-1e9
old_error=float(np.max(abs(p0[live].sum(axis=-1)-1)))
new_error=float(np.max(abs(p1[live].sum(axis=-1)-1)))
assert old_error>1e-12, 'Regression fixture did not reproduce cancellation'
assert new_error<1e-12 and np.isfinite(p1).all() and p1.min()>=0 and p1.max()<=1
np.testing.assert_array_equal(values0,values1)
np.testing.assert_array_equal(p1[~live],0.)
# Run the existing exact transformed solver tests with the candidate patch.
def proposed(text):
 transformed=original_patch(text)
 assert transformed.count(old)==1
 return transformed.replace(old,new,1)
patch_solver.patch_text=proposed
tests_solver.PATCHED=corrected
tests_solver.NS=tests_solver.make_namespace()
suite=unittest.defaultTestLoader.loadTestsFromTestCase(tests_solver.SolverPatchTests)
result=unittest.TextTestRunner(verbosity=2).run(suite)
assert result.wasSuccessful()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
receipt=dict(status='passed_proposed_numerical_repair_not_installed',model_solves=0,tests=result.testsRun,old_maximum_probability_sum_error=old_error,new_maximum_probability_sum_error=new_error,inclusive_values_exactly_unchanged=True,dead_menu_mask_unchanged=True,probability_tolerance_unchanged=1e-12,original_patch_sha256=sha(BASE/'patch_solver.py'),test_source_sha256=sha(__file__),candidate_source_sha256=sha(OUT/'solver_candidate.py'))
(OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
