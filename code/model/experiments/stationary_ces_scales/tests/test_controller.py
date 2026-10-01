"""Exercise the actual bounded optimizer loop without a lifecycle solve."""
import importlib.util,json,tempfile,time,unittest
from pathlib import Path
SRC=Path(__file__).resolve().parents[1]/'run_ces.py'
spec=importlib.util.spec_from_file_location('ces_controller',SRC);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Controller(unittest.TestCase):
    def test_exact_loop_zero_solves(self):
        seed={f'x{i}':.2 for i in range(10)};bounds={k:(0,1) for k in seed};seen=[]
        def evaluate(label,p,deadline):
            seen.append((label,p));return dict(status='passed',residual=[p[k]-.3 for k in seed],lifecycle_solves=0)
        with tempfile.TemporaryDirectory() as folder:
            r=m.controller(folder,seed,bounds,evaluate,evaluate,time.time()+3600,toy=True)
            self.assertEqual(r['lifecycle_solves'],0);self.assertLessEqual(r['completed_full_ge'],20)
            self.assertEqual(seen[0][1],seen[1][1]);self.assertEqual(seen[-1][1],r['selected']['parameters'])
            self.assertGreaterEqual(r['completed_full_ge'],13)
            self.assertEqual(r['status'],'selected_numerically_verified')
    def test_lambda_step(self):
        _,_,a,steps=m.simplex({'lambda_housing':.2},{'lambda_housing':(0,1)},('lambda_housing',))
        self.assertEqual(steps['lambda_housing'],.05);self.assertAlmostEqual(a[1,0],.25)
if __name__=='__main__':unittest.main()
