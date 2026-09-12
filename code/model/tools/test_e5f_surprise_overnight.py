import unittest
import numpy as np
from run_e5f_successive_surprises_overnight import next_psi
from run_e5f_successive_surprise_policy import policy_residual
class OvernightTest(unittest.TestCase):
 def test_secant_bracket_and_bounds(self):
  self.assertAlmostEqual(next_psi([(.1,-.2),(.2,.2)],.2,-.01,(0,.3)),.15)
  self.assertLessEqual(next_psi([(.1,.2),(.2,.3)],.2,-.01,(.05,.3)),.3)
  self.assertGreaterEqual(next_psi([(.1,.2),(.2,.3)],.2,-.01,(.05,.3)),.05)
 def test_fiscal_ledgers_are_separate_and_scaled(self):
  r=policy_residual(10,10,dict(payroll_tax_revenue=2,pension_outlays=2),1,.9,True)
  np.testing.assert_allclose(r,[0,0,20])
  r=policy_residual(10,10,dict(payroll_tax_revenue=2,pension_outlays=1.9),1,1,True)
  np.testing.assert_allclose(r,[0,10,0])
 def test_zero_unrebated_budget_has_no_transfer_root(self):
  self.assertEqual(policy_residual(1,1,dict(payroll_tax_revenue=0,pension_outlays=0),.4,0,False),[0,0])
 def test_three_coordinate_root_keeps_strict_physical_gates(self):
  import time
  from e5f_matched_pf_path_root import solve_price_path
  target=np.array([.7,2.,.3]);scale=np.array([1.,200.,200.])
  r=solve_price_path(initial_prices=target*.99,evaluate=lambda x:dict(residual=scale*np.log(target/x),mapping_valid=True),project=lambda x:x,
   slope=1.,market_tolerance=2e-4,max_log_step=.2,damping=1.,max_evaluations=8,deadline_monotonic=time.monotonic()+5,
   max_condition_number=1e10,worsening_factor=1.5,final_reproduction_tolerance=2e-10,default_jacobian=-np.diag(scale))
  self.assertTrue(r['converged']);np.testing.assert_allclose(r['final']['prices'],target,rtol=0,atol=1e-10)
if __name__=='__main__':unittest.main()
