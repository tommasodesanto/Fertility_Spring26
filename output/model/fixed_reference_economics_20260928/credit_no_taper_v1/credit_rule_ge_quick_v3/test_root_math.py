import math
import importlib.util
from pathlib import Path

spec=importlib.util.spec_from_file_location('ge',Path(__file__).with_name('run_credit_rule_ge_quick.py')); ge=importlib.util.module_from_spec(spec);spec.loader.exec_module(ge)

def test_closed_accounting_and_logsecant():
 parent=ge.load(ge.PARENT,'parent_for_test')
 c=parent.closed_accounting(2.,4.2,3.,12.,1.1)
 assert c['renewal_residual']==0 and c['population_scale']==4 and c['absolute_housing_residual']==0
 x=ge.logsecant({'factor':.8,'residual':-.2},{'factor':1.2,'residual':.2})
 assert .8 < x < 1.2 and math.isfinite(x)

def test_root_never_accepts_gate_failure():
 assert not ge.root_eligible({'residual':0.,'gate_passed':False})
 assert ge.root_eligible({'residual':0.,'gate_passed':True})
