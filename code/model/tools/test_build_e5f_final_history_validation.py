import importlib.util, json, tempfile, unittest, shutil
from pathlib import Path
P=Path(__file__).with_name('build_e5f_final_history_validation.py'); s=importlib.util.spec_from_file_location('a',P); m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
class T(unittest.TestCase):
 def setUp(self):
  self.d=Path(tempfile.mkdtemp()); self.p=self.d/'model_2023.json'; self.x={'calendar_year':2023,'fertility_stock_timing':{'parity_shares_40_44':{'1':.2,'2':.5,'3plus':.3},'moments':{'childless_rate_40_44':.1,'exactly_one_among_mothers_40_44':.22,'period_mean_age_first_birth':27.,'period_share_first_births_age30plus':.3}},'fertility':{'birth_flow_third_bin_entry':[1.],'birth_flow_explicit':[3.],'birth_flow_topcode_adjusted':[4.]},'housing_wealth':{'moments':{'aggregate_mean_occupied_rooms_capped9_18_85':3.,'own_rate_30_55':.6,'prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9':.4,'aggregate_wealth_to_annual_gross_labor_earnings':2.,'annual_bequest_flow_to_aggregate_wealth':.01,'old_total_wealth_to_annual_income_p90_p50_7684':3.}},'dated_first_birth_rooms':{'housing_response':.1},'recent_parent':{'model_value':.05}}; self.p.write_text(json.dumps(self.x)); self.e={'cps':{'tfr':2.,'childless_rate':.1,'parity_share_1':.2},'nchs':{'mean_age':27.,'share30':.3},'acs':{'mean_rooms':3.,'ownership':.6,'family_rooms':.4,'recent_parent':.05},'psid':{'first_birth_rooms':.1,'old_dispersion':3.},'wealth':2.,'bequest':.01}
 def test_complete_fixture(self): self.assertEqual(len(m.build_table(self.p,self.d/'o',self.e,{'fixture':'synthetic'},True)),13)
 def test_wrong_year(self): self.x['calendar_year']=2022; self.p.write_text(json.dumps(self.x)); self.assertRaises(ValueError,m.build_table,self.p,self.d/'o',self.e)
 def test_loader_metadata(self): self.assertIn('window == pooled_2005_2019',m.load_empirical()[1]['predicates']['wealth'])
 def tearDown(self): shutil.rmtree(self.d)
 def test_missing_failed_or_unlinked_verification(self):
  self.assertRaises(ValueError,m.build_table,self.p,self.d/'o')
  v=self.p.with_name('verification.json');v.write_text(json.dumps({'status':'FAIL'}))
  self.assertRaises(ValueError,m.build_table,self.p,self.d/'o')
  check={'status':'PASS','finite_converged':True,'verification_method':'native_saved_snapshot_aggregate_match','year':2023,'historical_fit_status':{'complete':True},'root_receipt_sha256':'a'*64}
  v.write_text(json.dumps(check));self.assertRaises(ValueError,m.build_table,self.p,self.d/'o')
  self.x['forecast_receipt_sha256']='a'*64;self.p.write_text(json.dumps(self.x))
  self.assertEqual(len(m.build_table(self.p,self.d/'o')),13)
 def test_fixture_is_labelled_and_no_hidden_override(self):
  self.assertRaises(ValueError,m.build_table,self.p,self.d/'o',self.e)
  m.build_table(self.p,self.d/'o',self.e,fixture=True)
  meta=json.loads((self.d/'o/validation_2023_manifest.json').read_text())
  self.assertTrue(meta['deterministic_fixture'])
if __name__=='__main__': unittest.main()
