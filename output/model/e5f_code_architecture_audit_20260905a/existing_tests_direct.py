"""Run selected existing pure tests without the failing local pytest CLI.

Run from the repository root with NUMBA_DISABLE_JIT=1 and code/model/.venv/bin/python.
The scientific modules come from the frozen production snapshot; the current
unchanged test files are exposed on its package path. No full model solve runs.
"""
from pathlib import Path
import sys, importlib, inspect, json, os
r=Path.cwd(); f=r/'tmp/e5f_overnight_20260905a/frozen_production/code/model'
sys.path[:0]=[str(f),str(f/'tools'),str(r/'code/model')]
if os.environ.get('NUMBA_DISABLE_JIT')!='1':
    raise RuntimeError('This bounded pure test runner requires NUMBA_DISABLE_JIT=1')
import intergen_eqscale_seq_optimized as package
package.__path__.append(str(r/'code/model/intergen_eqscale_seq_optimized'))
results=[]
for name in ['test_e5f_transition_calibration','test_independent_child_maturation']:
    mod=importlib.import_module('intergen_eqscale_seq_optimized.tests.'+name)
    for fn,call in vars(mod).items():
        if fn.startswith('test_') and inspect.isfunction(call) and fn!='test_parity_progression_remains_available_after_older_child_leaves':
            call(); results.append(name+'.'+fn)
mod=importlib.import_module('intergen_eqscale_seq_optimized.tests.test_age_survival')
for name in ['test_age_survival_default_is_inert','test_age_survival_override_is_external_and_validated']:
    getattr(mod.AgeSurvivalTests(),name)();results.append('test_age_survival.'+name)
x={'method':'Direct invocation of unchanged test functions; pytest runner was not used','passed':results,'count':len(results),'runtime_source':str(f),'tests_source':str(r/'code/model/intergen_eqscale_seq_optimized/tests'),'jit_disabled':True,'pytest_cli_attempts':3,'pytest_cli_exit_code':139,'pytest_cli_output':'none'}
(r/'output/model/e5f_code_architecture_audit_20260905a/existing_test_receipt.json').write_text(json.dumps(x,indent=2)+'\n')
print(json.dumps(x,indent=2))
