import json,os,shutil,tempfile,unittest
from pathlib import Path
from build_e5f_final_history_plots import build,sha,DEFAULT_DATA,csvread

class FinalHistoryPlotsTest(unittest.TestCase):
    def setUp(self):
        self.root=Path(tempfile.mkdtemp());self.case=self.root/'case';self.readout=self.root/'readout'
        self.case.mkdir();self.readout.mkdir();self.fits=[]
        (self.case/'contract_receipt.json').write_text(json.dumps(dict(case='A0',count=6)))
        for i,year in enumerate((2007,2011,2015,2019)):
            d=self.case/f'window_{year}/trial_00';d.mkdir(parents=True)
            receipt=dict(case='A0',count=6,start_year=year,psi=.1+i*.01,converged=True,
                finite_horizon_market_fiscal_converged=True,final={'mapping_valid':True},final_reproduction_max_abs=0.)
            (d/'root_receipt.json').write_text(json.dumps(receipt))
            rows=[dict(calendar_year=y,asset_price=1+i/10,renter_price=.2,housing_demand=5.,household_heads=1.,resident_persons=2.) for y in range(year,year+24,4)]
            fertility=[dict(calendar_year=y,period_tfr_topcode_adjusted=2-i*.1-j*.01) for j,y in enumerate(range(year,year+24,4))]
            (d/'rows.json').write_text(json.dumps(rows));(d/'fertility.json').write_text(json.dumps(fertility))
            self.fits.append(dict(year=year,model=2-i*.1,target=2-i*.1,gap=0.,psi=receipt['psi'],folder=str(d)))
        self.write_fit()
        ages=[int(r['age_lower']) for r in csvread(DEFAULT_DATA/'actual2023_age_housing_levels.csv')]
        profile=dict(totals=dict(households=1.,consumption=3.,rooms=5.),number_children_mass=[.2,.3,.3,.2],
            rows=[dict(age=a,households=.1,owners=.05,capped_rooms=.4,with_children=.02) for a in ages])
        last=self.case/'window_2019/trial_00/root_receipt.json';digest=sha(last)
        model=dict(calendar_year=2023,finite_converged=True,forecast_receipt_sha256=digest,profile=profile,deterministic_fixture=True,
            fertility=dict(age_cell_start=[18,22,26,30,34,38,42,46,50],age_specific_birth_rate_topcode_adjusted=[.01,.03,.08,.10,.07,.02,.01,.001,0.]))
        check=dict(status='PASS',finite_converged=True,verification_method='native_saved_snapshot_aggregate_match',
            historical_fit_status={'complete':True,'realized':self.fits},root_receipt_sha256=digest,horizon_verified=False)
        (self.readout/'model_2023.json').write_text(json.dumps(model));(self.readout/'verification.json').write_text(json.dumps(check))
    def write_fit(self):
        (self.case/'realized_fit.json').write_text(json.dumps(self.fits))
        (self.case/'finite_history_complete.json').write_text(json.dumps(dict(realized=self.fits,horizon_verified=False)))
    def tearDown(self):shutil.rmtree(self.root)
    def test_five_figures_and_correct_clock_and_data(self):
        out=self.root/'figures';m=build(self.case,self.readout,out=out)
        self.assertEqual(set(m['figures']),{'historical_fertility','prices_quantities_path','equilibrium_2023','lifecycle_2023','fertility_age_2023'})
        f=m['figures'];self.assertEqual(f['historical_fertility']['future_years'],[2023,2027,2031,2035,2039,2043])
        self.assertTrue(all(v==50 for v in f['lifecycle_2023']['series']['owners']['model']))
        self.assertTrue(f['lifecycle_2023']['series']['owners']['data'])
        self.assertTrue(f['fertility_age_2023']['data'])
        self.assertGreater((out/'e5f_final_history_figures.pdf').stat().st_size,1000)
        if os.environ.get('E5F_PLOT_FIXTURE_DIR'):
            target=Path(os.environ['E5F_PLOT_FIXTURE_DIR']);shutil.copytree(out,target,dirs_exist_ok=True)
    def test_bad_fit_and_readout_link_rejected(self):
        self.fits[0]['model']+=.01;self.fits[0]['gap']=.01;self.write_fit()
        with self.assertRaises(ValueError):build(self.case,self.readout)
        self.fits[0]['model']-=.01;self.fits[0]['gap']=0.;self.write_fit()
        model=json.loads((self.readout/'model_2023.json').read_text());model['forecast_receipt_sha256']='0'*64
        (self.readout/'model_2023.json').write_text(json.dumps(model))
        with self.assertRaises(ValueError):build(self.case,self.readout)
if __name__=='__main__':unittest.main()
