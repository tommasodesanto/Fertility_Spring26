import csv,hashlib,json,sys,tempfile,shutil
from pathlib import Path
from unittest import TestCase,main
sys.path.insert(0,str(Path(__file__).parent)); from build_e5f_final_policy_readout import build
def put(p,x): p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x))
def dg(p): return hashlib.sha256(p.read_bytes()).hexdigest()
class T(TestCase):
 def fx(self):
  z=Path(tempfile.mkdtemp());self.addCleanup(shutil.rmtree,z);c=z/'case';fit=[{'year':y,'psi':.4,'target':2.,'model':2.,'gap':0.} for y in (2007,2011,2015,2019)];put(c/'realized_fit.json',fit);put(c/'finite_history_complete.json',{'realized':fit});put(c/'contract_receipt.json',{'case':'A0_6','count':2})
  for n,t,o in [('baseline_rebate',.01,.5),('tax2_rebate',.02,.6)]:
   d=c/'policies'/n; rows=[dict(calendar_year=y,asset_price=1.,renter_price=.2,housing_demand=2.,resident_persons=1.,household_heads=1.,birth_children_topcode_adjusted=.2,owner_rate=o,pension_period_units=1.,equal_transfer_period_units=.1,relative_market_residual=0.,payroll_tax_revenue=1.,pension_outlays=1.,property_tax_revenue=1.,equal_transfer_outlays=1.) for y in (2023,2027)]
   put(d/'summary.json',{'finite_converged':True,'annual_tax':t});put(d/'root_receipt.json',{'case':'A0_6','count':2,'start_year':2023,'psi':.4,'converged':True,'finite_horizon_market_fiscal_converged':True,'horizon_verified':False,'production_eligible':False,'final_reproduction_max_abs':0.,'final':{'mapping_valid':True,'residual':[0.]*9}});put(d/'rows.json',rows);put(d/'fertility.json',[{'calendar_year':y,'period_tfr_topcode_adjusted':2.} for y in (2023,2027)])
  s=z/'proof.json';put(s,{'status':'passed','common_initial_g_pre':True,'common_grid':True,'common_supply_rule':True,'all_other_parameter_fields_exact':True,'worker_income_exact':True,'case_contract_sha256':dg(c/'contract_receipt.json'),'policy_root_sha256':{n:dg(c/'policies'/n/'root_receipt.json') for n in ('baseline_rebate','tax2_rebate')}});return c,s,z
 def test_build_percent_pp(self):
  c,s,z=self.fx();build(c,z/'o',s);r=next(x for x in csv.DictReader((z/'o/comparison.csv').read_text().splitlines()) if x['metric']=='owner_rate');self.assertAlmostEqual(float(r['absolute_change']),.1);self.assertAlmostEqual(float(r['percent_change']),20.);self.assertAlmostEqual(float(r['percentage_point_change']),10.)
 def bad(self,f):
  c,s,z=self.fx();f(c,s);self.assertRaises(ValueError,build,c,z/'o',s)
 def test_rejections(self):
  for k in ('year','psi','unconverged','replay','residual'):
   def f(c,s,k=k):
    p=c/'policies/tax2_rebate/root_receipt.json';x=json.loads(p.read_text());x['start_year']=2027 if k=='year' else x['start_year'];x['psi']=.3 if k=='psi' else x['psi'];x['converged']=False if k=='unconverged' else x['converged'];x['final_reproduction_max_abs']=3e-10 if k=='replay' else x['final_reproduction_max_abs'];x['final']['residual'][0]=3e-4 if k=='residual' else x['final']['residual'][0];put(p,x);q=json.loads(s.read_text());q['policy_root_sha256']['tax2_rebate']=dg(p);put(s,q)
   self.bad(f)
 def test_bad_state(self):
  c,s,z=self.fx();put(s,{'status':'failed','common_initial_g_pre':False});self.assertRaises(ValueError,build,c,z/'o',s)
 def test_changed_supply(self):
  c,s,z=self.fx();q=json.loads(s.read_text());q['common_supply_rule']=False;put(s,q);self.assertRaises(ValueError,build,c,z/'o',s)
 def test_conditional_history_label(self):
  c,s,z=self.fx();p=c/'contract_receipt.json';q=json.loads(p.read_text());q.update(conditional_history_count=6,history_refitted=False);put(p,q)
  (c/'finite_history_complete.json').rename(c/'conditioning_history_complete.json')
  q=json.loads(s.read_text());q.update(case_contract_sha256=dg(p),conditional_history_count=6);put(s,q)
  receipt=build(c,z/'o',s);self.assertEqual(receipt['conditional_history_count'],6)
  q['conditional_history_count']=24;put(s,q);self.assertRaises(ValueError,build,c,z/'bad',s)
if __name__=='__main__':main()
