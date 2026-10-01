"""Bounded pure checks; no native model evaluations."""
import json,tempfile,time,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import inputs,runner
class Checks(unittest.TestCase):
 def test_18_starts(self):
  self.assertEqual(len(inputs.LANES),18)
  for lane,cfg in inputs.LANES.items():
   inputs.check_point(cfg['seed'],cfg['bounds']);self.assertEqual(len(cfg['seed']),8 if cfg['arm']=='constant_alpha' else 9)
   self.assertEqual(cfg['seed_provenance']['full_ge_repeats'],2)
  self.assertEqual(inputs.LANES['floor_s0']['seed']['h_P'],1.8900476600128304)
  for a in ('floor','no_A','constant_alpha'):
   self.assertEqual(inputs.LANES[a+'_s4']['seed']['first_birth_fixed_cost'],0.)
   self.assertEqual(inputs.LANES[a+'_s5']['seed']['first_birth_fixed_cost'],1.6023128423441799)
 def test_entry_and_binding(self):
  for arm in ('floor','no_A','constant_alpha'):
   lane=arm+'_s0';P,g=inputs.proposal(lane);P,r=inputs.entry(P,g,'nonnegative_mean')
   self.assertEqual(r['negative_wealth_share'],0.);self.assertEqual(r['clipped_draw_mass'],0.)
   seed,bounds,_=inputs.seed_and_bounds(lane);Q=inputs.bind(P,seed,bounds,arm);expected=runner.expected_parameters(seed,arm=arm)
   self.assertEqual(Q.hbar_first_child_jump,expected['h_P']);self.assertEqual(Q.delta_alpha_jump,expected['delta_alpha_jump']);self.assertEqual(Q.psi_child,P.psi_child)
   np.testing.assert_array_equal(Q.H0,P.H0)
   with self.assertRaisesRegex(RuntimeError,'must be nonnegative_mean'):inputs.entry(P,g,arm)
 def test_search_actual_loop_8_and_9(self):
  for arm in ('floor','no_A','constant_alpha'):
   lane=arm+'_s0';seed,bounds,_=inputs.seed_and_bounds(lane)
   with tempfile.TemporaryDirectory() as d:
    out=Path(d);r=runner.search(out,seed,bounds,runner.mocked_evaluator(out,seed,bounds,lane),time.time()+5400,mock=True,lane=lane)
    self.assertEqual(r['completed_full_ge'],32);self.assertEqual(r['selected_repeats'],2);self.assertEqual(r['lifecycle_solves'],0)
    rounds=json.loads((out/'rounds.json').read_text());self.assertTrue(all(x['rank']==len(seed) for x in rounds[:-1]));self.assertEqual(rounds[-1]['status'],'incomplete_jacobian')
 def test_failed_baseline_fatal(self):
  seed,bounds,_=inputs.seed_and_bounds('floor_s0')
  with tempfile.TemporaryDirectory() as d:
   with self.assertRaisesRegex(RuntimeError,'Baseline full GE failed'):runner.search(Path(d),seed,bounds,lambda *_:dict(status='inadmissible_numerical'),time.time()+5400,mock=True)
 def test_failed_probe_not_zero_derivative(self):
  for arm in ('floor','constant_alpha'):
   lane=arm+'_s0';seed,bounds,_=inputs.seed_and_bounds(lane)
   with tempfile.TemporaryDirectory() as d:
    out=Path(d);good=runner.mocked_evaluator(out,seed,bounds,lane)
    def ev(label,p,end):return dict(status='inadmissible_numerical',lifecycle_solves=0) if 'probe' in label else good(label,p,end)
    r=runner.search(out,seed,bounds,ev,time.time()+5400,mock=True,lane=lane)
    self.assertEqual(r['search_stop_reason'],'incomplete_jacobian');self.assertEqual(r['selected_repeats'],2);self.assertNotIn('rank',r['identification'])
 def test_rank_deficiency(self):
  for arm in ('floor','constant_alpha'):
   lane=arm+'_s0';seed,bounds,_=inputs.seed_and_bounds(lane)
   with tempfile.TemporaryDirectory() as d:
    r=runner.search(Path(d),seed,bounds,lambda *_:dict(status='passed',residual=[1.]*10,lifecycle_solves=0),time.time()+5400,mock=True,lane=lane)
    self.assertEqual(r['search_stop_reason'],'rank_deficient');self.assertEqual(r['completed_full_ge'],len(seed)+4)
 def test_reserve_and_budget_exit(self):
  self.assertEqual(runner.reserve_seconds('floor_s0',210),1256.)
  self.assertGreaterEqual(runner.reserve_seconds('floor_s0',800),2790.)
  seed,bounds,_=inputs.seed_and_bounds('floor_s0')
  with tempfile.TemporaryDirectory() as d:
   out=Path(d);good=runner.mocked_evaluator(out,seed,bounds,'floor_s0');ends=[];deadline=time.time()+5400
   def ev(label,p,end):
    ends.append((label,end));return dict(status='budget_exhausted',lifecycle_solves=0) if 'probe' in label else good(label,p,end)
   with patch.object(runner,'compare_repeated',return_value={}):r=runner.search(out,seed,bounds,ev,deadline,mock=False)
   self.assertEqual(r['selected_repeats'],2);self.assertTrue(all(e==deadline for k,e in ends if 'selected_repeat' in k))
 def test_actual_target_and_all31_expected(self):
  self.assertEqual(len(runner.PLAN['target_contract']),14);self.assertEqual(sum(x['role']=='scored' for x in runner.PLAN['target_contract']),10)
  for lane,cfg in inputs.LANES.items():
   expected=runner.expected_parameters(cfg['seed'],arm=cfg['arm']);self.assertEqual(len(expected),31)
   self.assertEqual(expected['delta_alpha'],0.);self.assertEqual(expected['delta_alpha_jump'],cfg['seed'].get('delta_alpha_jump',0.))
if __name__=='__main__':unittest.main()
