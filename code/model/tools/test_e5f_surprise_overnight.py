import unittest
import numpy as np
from run_e5f_successive_surprises_overnight import next_psi,attempt_forecast,restore_fitted_prefix,sha
from run_e5f_successive_surprise_policy import policy_residual
class OvernightTest(unittest.TestCase):
 def test_rejected_numerical_trial_does_not_end_candidate_loop(self):
  import tempfile,json
  from pathlib import Path
  with tempfile.TemporaryDirectory() as folder:
   p=Path(folder)/'rejected.json';accepted=[]
   def bad():raise RuntimeError('population mass gate')
   for solve in (bad,lambda:42):
    result=attempt_forecast(solve,p)
    if result is not None:accepted.append(result)
   self.assertEqual(accepted,[42]);self.assertEqual(json.loads(p.read_text())['error'],'population mass gate')
   def invalid():raise ValueError('wrong contract')
   with self.assertRaises(ValueError):attempt_forecast(invalid,p)
 def test_resume_requires_pinned_fitted_prefix_and_matching_clock(self):
  import tempfile,json,gzip,pickle,copy
  from pathlib import Path
  from types import SimpleNamespace
  with tempfile.TemporaryDirectory() as folder:
   root=Path(folder)
   plan=dict(terminal_template={},source_root='fixture',target_fingerprint='fixed',history_root_controls={},
    outside_origin_entry_share=.169,fertility_fit_tolerance=.005,forecast_counts=[6],finite_sequence_diagnostic=True,skip_policies=True)
   values=dict(source_plan=copy.deepcopy(plan),realized_fit=[dict(year=2007,error_abs=.001,psi=.14)],
    receipt=dict(finite_horizon_market_fiscal_converged=True,start_year=2007,psi=.14))
   for key,value in values.items():(root/key).write_text(json.dumps(value))
   with gzip.open(root/'checkpoint','wb') as f:pickle.dump(SimpleNamespace(year=2011),f)
   plan['resume_fitted_prefix']={key:dict(path=str(root/key),sha256=sha(root/key)) for key in (*values,'checkpoint')}
   inherited,rows,receipt=restore_fitted_prefix(plan)
   self.assertEqual(inherited.year,2011);self.assertEqual(len(rows),1)
   altered=copy.deepcopy(plan);altered['history_root_controls']={'fiscal_tolerance':.01}
   with self.assertRaises(ValueError):restore_fitted_prefix(altered)
   (root/'realized_fit').write_text('[]')
   with self.assertRaises(ValueError):restore_fitted_prefix(plan)
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
