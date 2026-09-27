import math, unittest, tempfile, csv
from pathlib import Path
from run_e5f_local_jacobian import collect, write
from run_e5f_local_jacobian import make_points
class TestPoints(unittest.TestCase):
 def test_central_geometry(self):
  p={'H0':6.,'beta_annual':.96,'delta_alpha_jump':.14,'child_benefit_curvature':.07}
  rs=[dict(parameter=k,lower=0,upper=100) for k in p]
  rows=make_points(p,rs)
  self.assertEqual(len(rows),8)
  self.assertAlmostEqual(math.log(rows[1]['point']['H0']/rows[0]['point']['H0']),.04)
  self.assertAlmostEqual(rows[3]['point']['beta_annual']-rows[2]['point']['beta_annual'],.001)
  self.assertAlmostEqual(rows[7]['point']['child_benefit_curvature']-rows[6]['point']['child_benefit_curvature'],.016)
 def test_reject_bound_crossing(self):
  with self.assertRaises(AssertionError):make_points({'H0':1},[dict(parameter='H0',lower=1,upper=2)])
 def test_known_weighted_derivative(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);anchor=root/'anchor';anchor.mkdir()
   fields=['moment','target','model','gap','weight','loss_contribution']
   def fit(folder,model):
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/'target_fit.csv').open('w',newline='') as f:
     w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerow(dict(moment='test',target=1,model=model,gap=model-1,weight=4,loss_contribution=4*(model-1)**2))
   fit(anchor,2);fit(root/'a_minus/case',1);fit(root/'a_plus/case',3)
   write(root/'a_minus/case/receipt.json',{'normalization':{'psi_child':1}})
   write(root/'a_plus/case/receipt.json',{'normalization':{'psi_child':2}})
   plan=dict(anchor_case=str(anchor),point={'a':2},cases=[dict(id='a_minus',parameter='a',h=.02,coordinate='log'),dict(id='a_plus',parameter='a',h=.02,coordinate='log')])
   result=collect(plan,root)
   self.assertAlmostEqual(result['raw_jacobian'][0][0],50)
   self.assertAlmostEqual(result['weighted_jacobian'][0][0],100)
   self.assertAlmostEqual(result['singular_values'][0],100)
   self.assertAlmostEqual(result['loss_gradient'][0],400)
   self.assertAlmostEqual(result['normalized_benefit_derivatives'][0],25)
if __name__=='__main__':unittest.main()
