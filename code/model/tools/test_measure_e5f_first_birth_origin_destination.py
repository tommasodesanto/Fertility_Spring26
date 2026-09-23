"""Small exact-loop tests for the read-only active birth observer."""
from types import SimpleNamespace
import unittest
import numpy as np
from measure_e5f_first_birth_origin_destination import assert_profile, origin_branches, verified_gap

class Native:
    @staticmethod
    def get_fecundity_by_age(P): return np.array([.5, 0.])
    @staticmethod
    def readiness_settled_state(P): return 0

class TestObserver(unittest.TestCase):
    def P(self, **kw):
        d=dict(housing_event_horizon=1,period_years=4,da=4,sequential_births=True,joint_nested_choice=False,fertility_nest_choice=False,two_shock_choice=False,child_state_mode='independent_count',J=2,A_f_start=1,A_f_end=1)
        d.update(kw); return SimpleNamespace(**d)
    def test_exact_first_birth_loop_and_one_block_smoke(self):
        g=np.zeros((1,1,1,2,2,4,4)); g[0,0,0,0,:,0,0]=[4.,8.]
        fp=np.zeros((1,1,1,2,2,4)); fp[0,0,0,0,:,1]=[.25,.5]
        e=SimpleNamespace(g_pre=g,policy=SimpleNamespace(fert_probs=fp)); P=self.P()
        tr,co,settled,counts=origin_branches(e,P,Native(),False)
        self.assertEqual(settled,0); self.assertAlmostEqual(tr.sum(),2.5); self.assertAlmostEqual(co.sum(),2.5); self.assertEqual(counts['positive_age_income_blocks'],2)
        tr,co,_,counts=origin_branches(e,P,Native(),True)
        self.assertAlmostEqual(tr.sum(),2.); self.assertAlmostEqual(co.sum(),2.); self.assertEqual(counts['selected_blocks'],1)
    def test_joint_profile_fails_closed(self):
        with self.assertRaisesRegex(ValueError,'without a joint nest'): assert_profile(self.P(joint_nested_choice=True))
    def test_reproduction_gate_rejects_mismatch_and_nonfinite_values(self):
        self.assertEqual(verified_gap(1.,1.,2e-10), 0.)
        with self.assertRaises(RuntimeError): verified_gap(1.,1.1,2e-10)
        for observed, expected, tolerance in [(float('nan'),1.,2e-10),(1.,float('inf'),2e-10),(1.,1.,float('nan')),(1.,1.,0)]:
            with self.assertRaises(ValueError): verified_gap(observed, expected, tolerance)
if __name__=='__main__': unittest.main()
