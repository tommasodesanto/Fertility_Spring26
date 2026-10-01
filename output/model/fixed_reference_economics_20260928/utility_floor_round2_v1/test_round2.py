"""Pure exact controller checks for the four-hour floor continuation."""
import json,tempfile,time,unittest
from pathlib import Path
from unittest.mock import patch
import inputs,runner
class Checks(unittest.TestCase):
 def test_eight_bounded_starts(self):
  self.assertEqual(len(inputs.LANES),8)
  for c in inputs.LANES.values():inputs.check_point(c['seed'],c['bounds']);self.assertEqual(len(c['seed']),9);self.assertEqual(c['arm'],'floor')
 def test_unchanged_targets_and_bounds(self):
  old=json.loads((inputs.HERE.parent/'utility_calibration_round1_v1/plan.json').read_text())
  self.assertEqual(runner.PLAN['target_contract'],old['target_contract'])
  self.assertEqual(runner.PLAN['lanes']['floor_s0']['bounds'],old['lanes']['floor_s0']['bounds'])
 def test_exact80_GE_loops(self):
  for lane in ('floor_s0','floor_s1'):
   seed,bounds,_=inputs.seed_and_bounds(lane)
   with tempfile.TemporaryDirectory() as d:
    out=Path(d);r=runner.search(out,seed,bounds,runner.mocked_evaluator(out,seed,bounds,lane),time.time()+14400,mock=True,lane=lane)
    self.assertEqual(r['completed_full_ge'],80);self.assertEqual(r['selected_repeats'],2);self.assertEqual(r['lifecycle_solves'],0)
 def test_incomplete_probe_and_fatal_baseline(self):
  seed,bounds,_=inputs.seed_and_bounds('floor_s0')
  with tempfile.TemporaryDirectory() as d:
   with self.assertRaisesRegex(RuntimeError,'Baseline'):runner.search(Path(d),seed,bounds,lambda *_:dict(status='inadmissible_numerical'),time.time()+14400,mock=True)
  with tempfile.TemporaryDirectory() as d:
   out=Path(d);good=runner.mocked_evaluator(out,seed,bounds,'floor_s0')
   def ev(label,p,end):return dict(status='inadmissible_numerical',lifecycle_solves=0) if 'probe' in label else good(label,p,end)
   r=runner.search(out,seed,bounds,ev,time.time()+14400,mock=True)
   self.assertEqual(r['search_stop_reason'],'incomplete_jacobian');self.assertNotIn('rank',r['identification']);self.assertEqual(r['selected_repeats'],2)
 def test_rank_and_accounting_failure(self):
  seed,bounds,_=inputs.seed_and_bounds('floor_s0')
  with tempfile.TemporaryDirectory() as d:
   r=runner.search(Path(d),seed,bounds,lambda *_:dict(status='passed',residual=[1.]*10,lifecycle_solves=0),time.time()+14400,mock=True)
   self.assertEqual(r['search_stop_reason'],'rank_deficient')
  with tempfile.TemporaryDirectory() as d:
   def ev(*_):raise RuntimeError('Accounting drift')
   with self.assertRaisesRegex(RuntimeError,'Accounting drift'):runner.search(Path(d),seed,bounds,ev,time.time()+14400,mock=True)
if __name__=='__main__':unittest.main()
