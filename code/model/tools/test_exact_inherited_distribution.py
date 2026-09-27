"""Pure gate fixtures: no model import or equilibrium solve."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
import numpy as np


def load_gate(model):
    source=Path(__file__).with_name('run_dynamic_population_transition.py').read_text()
    names={'InheritedDistributionInfeasible','_require_exact_inherited_distribution','gate_pre_fertility_distribution'}
    nodes=[n for n in ast.parse(source).body if getattr(n,'name',None) in names]
    future=ast.ImportFrom(module='__future__',names=[ast.alias(name='annotations')],level=0)
    module=ast.fix_missing_locations(ast.Module(body=[future]+nodes,type_ignores=[]))
    ns=dict(np=np,hashlib=hashlib,json=json,Path=Path,time=time,model=model)
    exec(compile(module,'<pure inherited gate>','exec'),ns)
    return ns


class ExactInheritedTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.ns=load_gate(SimpleNamespace(DEAD_MASS_TOL=1e-12))
        self.P=SimpleNamespace(J=1,age_start=18,da=4,user_cost_rate=.1,
            native_exact_inherited_distribution=True,native_inherited_distribution_evidence_dir=self.tmp.name)
        self.g=np.zeros((2,1,1,1,1,1,1));self.g[0]=.3;self.g[1]=.7
        self.policy=SimpleNamespace(V=np.zeros_like(self.g),price=np.array([1.]))
        self.grid=np.array([-1.,1.])

    def tearDown(self):self.tmp.cleanup()
    def gate(self):return self.ns['gate_pre_fertility_distribution'](self.g,self.policy,self.P,self.grid,None)

    def test_exact_array_preserved_without_aliasing(self):
        before=self.g.tobytes();out,projection=self.gate()
        self.assertEqual(out.tobytes(),before);self.assertEqual(projection,0.)
        self.assertIsNot(out,self.g);self.assertEqual(self.g.tobytes(),before)

    def test_above_existing_tolerance_rejected_and_recorded(self):
        self.g[0]=1e-8;self.policy.V[0]=-1e9
        before=self.g.tobytes()
        with self.assertRaises(self.ns['InheritedDistributionInfeasible']) as caught:self.gate()
        exc=caught.exception
        self.assertEqual(exc.classification,'inherited_distribution_infeasible')
        self.assertEqual(exc.dead_mass,1e-8)
        self.assertEqual(self.g.tobytes(),before)
        saved=json.loads(Path(exc.evidence_path).read_text())
        self.assertEqual(saved['census'][0]['wealth'],-1.)
        self.assertEqual(saved['census'][0]['value'],-1e9)
        self.assertEqual(saved['origin_sha256'],hashlib.sha256(self.g.tobytes()).hexdigest())
        self.assertFalse(saved['distribution_modified'])
        self.assertIsInstance(exc,RuntimeError)
        self.assertNotIn('exceeds',str(exc))

    def test_small_tail_preserved_and_reported_without_redistribution(self):
        self.g[0]=1e-14;self.policy.V[0]=-1e9
        before=self.g.tobytes();out,projection=self.gate()
        self.assertEqual(out.tobytes(),before);self.assertEqual(projection,0.)
        saved=json.loads(next(Path(self.tmp.name).glob('*.json')).read_text())
        self.assertEqual(saved['status'],'retained_below_existing_feasibility_tolerance')
        self.assertEqual(saved['feasibility_mass_tolerance'],1e-12)
        self.assertFalse(saved['distribution_modified'])

    def test_zero_mass_dead_cell_allowed(self):
        self.g[0]=0.;self.policy.V[0]=-1e10
        out,projection=self.gate();np.testing.assert_array_equal(out,self.g)
        self.assertEqual(projection,0.)

    def test_occupied_nonfinite_value_rejected(self):
        self.policy.V[0]=np.nan
        with self.assertRaises(self.ns['InheritedDistributionInfeasible']) as caught:self.gate()
        self.assertIsNone(caught.exception.census[0]['value'])

    def test_legacy_projection_and_native_gate_still_called(self):
        calls=[]
        def censor(g,v):
            calls.append('censor');mass=float(g[0].sum());g[1]+=g[0];g[0]=0;return mass
        def gate(*args,**kwargs):calls.append('native_gate')
        ns=load_gate(SimpleNamespace(_censor_entry_dead_mass=censor,_gate_dead_mass_at_age=gate))
        del self.P.native_exact_inherited_distribution
        out,projection=ns['gate_pre_fertility_distribution'](self.g,self.policy,self.P,self.grid,None)
        self.assertEqual(calls,['censor','native_gate']);self.assertEqual(projection,.3)
        self.assertEqual(float(out[1].sum()),1.);self.assertEqual(float(self.g[0].sum()),.3)

    def test_invalid_distribution_is_not_solver_infeasibility(self):
        self.g[0]=-1e-30
        with self.assertRaises(ValueError):self.gate()

if __name__=='__main__':unittest.main()
