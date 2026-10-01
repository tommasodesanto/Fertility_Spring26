import importlib.util
from pathlib import Path
import unittest
import numpy as np
s=importlib.util.spec_from_file_location('fast_controller',Path(__file__).with_name('controller.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Tests(unittest.TestCase):
 def test_original_formula_and_three_dampings(self):
  coordinates=tuple('abcdefghi');seed={k:.4 for k in coordinates};bounds={k:(0.,1.) for k in coordinates};steps={k:.01 for k in coordinates};rr=np.arange(10.)/10
  base=dict(label='000_baseline',parameters=seed,residual=rr.tolist());J=np.vstack((np.eye(9),np.ones(9)))
  probes=[dict(residual=(rr+.01*J[:,i]).tolist()) for i in range(9)]
  receipt,points=m.jacobian(None,base,probes,coordinates,steps,bounds)
  sv=np.linalg.svd(J,compute_uv=False);ridge=max(float(sv[0]**2)*1e-4,1e-10);delta=np.clip(np.linalg.solve(J.T@J+ridge*np.eye(9),-J.T@rr),-.03,.03)
  self.assertEqual(receipt['rank'],9)
  np.testing.assert_allclose(receipt['normalized_step'],delta,rtol=0,atol=1e-15)
  for p,d in zip(points,(.5,.2,1.)):np.testing.assert_allclose([p[k] for k in coordinates],.4+d*delta,rtol=0,atol=1e-15)
  self.assertEqual(receipt['center_parameters'],seed)
  self.assertFalse(receipt['fresh_at_selected'])
 def test_rank_deficient_blocks_GN(self):
  keys=tuple('abcdefghi');base=dict(label='base',parameters={k:.4 for k in keys},residual=[0.]*10)
  receipt,points=m.jacobian(None,base,[dict(residual=[0.]*10) for _ in keys],keys,{k:.01 for k in keys},{k:(0.,1.) for k in keys})
  self.assertEqual(receipt['rank'],0);self.assertEqual(points,[])
if __name__=='__main__':unittest.main()
