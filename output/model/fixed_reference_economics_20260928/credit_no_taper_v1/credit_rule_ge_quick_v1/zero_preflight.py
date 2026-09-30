"""No-lifecycle mechanics preflight for the staged credit-rule packet."""
import importlib.util
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod

ge=load('credit_ge_quick',HERE/'run_credit_rule_ge_quick.py')
api=load('credit_rule_api',HERE/'rule_interface.py')
class P:
    J=3; lambda_d=1.; survival_probs=np.array([1., .9, 1.]); use_age_survival=True
class M: pass
p=P(); api.apply_rule(M(),p,'author'); assert np.all(p.debt_taper_weights==0) and np.all(p.debt_caps==0)
p=P(); api.apply_rule(M(),p,'ours'); assert p.debt_taper_weights.tolist()==[0.,1.,0.,0.]
assert not ge.root_eligible({'residual':0.,'gate_passed':False})
assert ge.root_eligible({'residual':0.,'gate_passed':True})
print('zero preflight passed')
