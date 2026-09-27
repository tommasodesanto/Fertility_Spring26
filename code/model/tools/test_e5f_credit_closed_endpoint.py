import unittest
import math
from types import SimpleNamespace
import numpy as np
import run_e5f_credit_closed_endpoint as d

class EndpointTests(unittest.TestCase):
    def test_reference_closure_identity(self):
        r=d.closed_accounting(.1,.21,5.,10.,2.)
        self.assertAlmostEqual(r['renewal_residual'],0.)
        self.assertEqual(r['population_scale'],2.)
        self.assertEqual(r['absolute_housing_demand'],10.)
        self.assertAlmostEqual(r['absolute_adjusted_births']/2.1,r['absolute_entry'])

    def test_scale_does_not_repair_birth_imbalance(self):
        a=d.closed_accounting(.1,.1,5.,10.,2.)
        b=d.closed_accounting(.1,.1,5.,100.,2.)
        self.assertEqual(a['renewal_residual'],b['renewal_residual'])
        self.assertNotEqual(a['population_scale'],b['population_scale'])

    def test_root_budget_and_fixed_parameter_callback(self):
        calls=[]
        def evaluate(p):
            calls.append(p);return dict(price=p,renewal_residual=math.log(p/1.3))
        result=d.bounded_root(evaluate,1.)
        self.assertEqual(result['status'],'root_tolerance_met')
        self.assertLessEqual(len(calls),16)
        self.assertLess(abs(result['selected']['renewal_residual']),d.ROOT_TOL)

    def test_unbracketed_does_not_claim_root(self):
        result=d.bounded_root(lambda p:dict(price=p,renewal_residual=1.),1.)
        self.assertEqual(result['status'],'no_sign_change_on_declared_schedule')
        self.assertEqual(len(result['rows']),5)

    def test_multiple_roots_not_selected(self):
        result=d.bounded_root(lambda p:dict(price=p,renewal_residual=(p-.8)*(p-1.2)),1.)
        self.assertEqual(result['status'],'multiple_candidate_roots_not_selected')
        self.assertIsNone(result['selected'])

    def test_invalid_accounting(self):
        for values in [(0,1,1,1,1),(.1,.2,0,1,1),(.1,.2,1,-1,1)]:
            with self.assertRaises(ValueError):d.closed_accounting(*values)

    def test_parameters_snapshot_no_mutation(self):
        P=SimpleNamespace(q=.08243216,period_years=4,phi=np.array([.8]),psi_child=.144)
        snapshot=d.primitive_snapshot(P)
        self.assertAlmostEqual(snapshot['annual_q'],.02)
        self.assertEqual(snapshot['phi'],[.8])
        self.assertEqual(P.psi_child,.144)

if __name__=='__main__':unittest.main()
